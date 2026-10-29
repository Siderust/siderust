// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-kernel support.
//!
//! ## Scientific scope
//!
//! This module parses frame metadata needed to map NAIF frame identifiers to
//! inertial, PCK-backed, CK-backed, and text-kernel frame specifications.
//!
//! ## Technical scope
//!
//! Text FK kernels support metadata lookup and fixed MATRIX/ANGLES TK
//! frames. Dynamic frame specifications are not evaluated.
//!
//! ## References
//!
//! - NAIF. *Frames Required Reading*.
//! - NAIF. *Kernel Required Reading*.

use super::text::{TextKernel, TextValue};
use super::SpiceError;
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar};
use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

/// SPICE frame class code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameClass {
    /// Built-in or FK-declared inertial frame.
    Inertial = 1,
    /// Body-fixed frame whose orientation is driven by PCK data.
    Pck = 2,
    /// Spacecraft/instrument frame whose orientation is driven by CK data.
    Ck = 3,
    /// Fixed text-kernel frame class code 4 (legacy variant name).
    TwoVector = 4,
    /// Dynamic frame evaluated from a runtime definition.
    Dynamic = 5,
}

/// TK-frame specification.
#[derive(Debug, Clone)]
pub enum TkSpec {
    /// Constant 3×3 rotation matrix, row-major, from defined to relative frame.
    Matrix([[f64; 3]; 3]),
}

/// One frame definition parsed from an FK kernel.
#[derive(Debug, Clone)]
pub struct FrameSpec {
    /// NAIF integer frame ID.
    pub id: i32,
    /// Frame name (uppercase).
    pub name: String,
    /// Frame class.
    pub class: FrameClass,
    /// Class-specific ID (for example a body or instrument ID).
    pub class_id: i32,
    /// Center body NAIF ID.
    pub center_id: i32,
    /// For TK frames: the constant rotation matrix and base frame name.
    pub tk_spec: Option<(TkSpec, String)>,
}

/// Parsed frame-kernel contents.
#[derive(Debug, Clone, Default)]
pub struct FrameKernel {
    /// All parsed frame definitions.
    pub frames: Vec<FrameSpec>,
}

impl FrameKernel {
    /// Parse an FK text kernel.
    pub fn from_text(src: &str) -> Result<Self, SpiceError> {
        let kernel = TextKernel::parse(src)?;
        let mut ids = Vec::new();
        for key in kernel.data.keys() {
            if let Some(id) = key
                .strip_prefix("FRAME_")
                .and_then(|rest| rest.strip_suffix("_NAME"))
                .and_then(|rest| rest.parse::<i32>().ok())
            {
                ids.push(id);
            }
        }
        ids.sort_unstable();
        ids.dedup();

        let mut frames = Vec::with_capacity(ids.len());
        for id in ids {
            let name = text_value(
                kernel.get(&format!("FRAME_{id}_NAME")).ok_or_else(|| {
                    SpiceError::FormatParse(format!("FK missing FRAME_{id}_NAME"))
                })?,
                &format!("FRAME_{id}_NAME"),
            )?
            .to_ascii_uppercase();
            let class = parse_frame_class(int_value(
                kernel.get(&format!("FRAME_{id}_CLASS")).ok_or_else(|| {
                    SpiceError::FormatParse(format!("FK missing FRAME_{id}_CLASS"))
                })?,
                &format!("FRAME_{id}_CLASS"),
            )?)?;
            let class_id = int_value(
                kernel.get(&format!("FRAME_{id}_CLASS_ID")).ok_or_else(|| {
                    SpiceError::FormatParse(format!("FK missing FRAME_{id}_CLASS_ID"))
                })?,
                &format!("FRAME_{id}_CLASS_ID"),
            )?;
            let center_id = int_value(
                kernel.get(&format!("FRAME_{id}_CENTER")).ok_or_else(|| {
                    SpiceError::FormatParse(format!("FK missing FRAME_{id}_CENTER"))
                })?,
                &format!("FRAME_{id}_CENTER"),
            )?;

            let tk_spec = if class == FrameClass::TwoVector {
                let spec_key = format!("TKFRAME_{id}_SPEC");
                let spec = kernel
                    .get(&spec_key)
                    .map(|value| text_value(value, &spec_key))
                    .transpose()?;
                match spec {
                    Some("MATRIX") => {
                        let matrix_key = format!("TKFRAME_{id}_MATRIX");
                        let relative_key = format!("TKFRAME_{id}_RELATIVE");
                        let values = kernel.get_f64_array(&matrix_key).ok_or_else(|| {
                            SpiceError::FormatParse(format!("FK missing or invalid {matrix_key}"))
                        })?;
                        if values.len() != 9 || values.iter().any(|v| !v.is_finite()) {
                            return Err(SpiceError::FormatParse(format!(
                                "{matrix_key} must contain 9 values"
                            )));
                        }
                        let relative = text_value(
                            kernel.get(&relative_key).ok_or_else(|| {
                                SpiceError::FormatParse(format!("FK missing {relative_key}"))
                            })?,
                            &relative_key,
                        )?
                        .to_ascii_uppercase();
                        Some((
                            TkSpec::Matrix([
                                [values[0], values[1], values[2]],
                                [values[3], values[4], values[5]],
                                [values[6], values[7], values[8]],
                            ]),
                            relative,
                        ))
                    }
                    Some("ANGLES") => {
                        let key = |suffix: &str| format!("TKFRAME_{id}_{suffix}");
                        let angles = kernel
                            .get_f64_array(&key("ANGLES"))
                            .ok_or_else(|| SpiceError::FormatParse("missing FK ANGLES".into()))?;
                        let axes = kernel
                            .get_f64_array(&key("AXES"))
                            .ok_or_else(|| SpiceError::FormatParse("missing FK AXES".into()))?;
                        let units = text_value(
                            kernel.get(&key("UNITS")).ok_or_else(|| {
                                SpiceError::FormatParse("missing FK UNITS".into())
                            })?,
                            &key("UNITS"),
                        )?;
                        let factor = match units {
                            "RADIANS" => 1.0,
                            "DEGREES" => core::f64::consts::PI / 180.0,
                            "ARCSECONDS" => core::f64::consts::PI / 648000.0,
                            _ => {
                                return Err(SpiceError::FormatParse(format!(
                                    "unsupported FK angle units {units}"
                                )))
                            }
                        };
                        if angles.len() != 3 || axes.len() != 3 {
                            return Err(SpiceError::FormatParse(
                                "FK ANGLES/AXES require three entries".into(),
                            ));
                        }
                        let mut matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
                        for (angle, axis) in angles.iter().zip(axes.iter()) {
                            if !angle.is_finite() || ![1.0, 2.0, 3.0].contains(axis) {
                                return Err(SpiceError::FormatParse(
                                    "invalid FK angle or axis".into(),
                                ));
                            }
                            matrix = super::bpc::mul(
                                matrix,
                                super::bpc::axis_rotation(*axis as i32, angle * factor),
                            );
                        }
                        let relative = text_value(
                            kernel.get(&key("RELATIVE")).ok_or_else(|| {
                                SpiceError::FormatParse("missing FK RELATIVE".into())
                            })?,
                            &key("RELATIVE"),
                        )?
                        .to_ascii_uppercase();
                        Some((TkSpec::Matrix(matrix), relative))
                    }
                    Some(other) => {
                        return Err(SpiceError::FormatParse(format!(
                            "unsupported TKFRAME spec '{other}' for frame {id}"
                        )));
                    }
                    None => None,
                }
            } else {
                None
            };

            if let Some((TkSpec::Matrix(m), _)) = &tk_spec {
                for i in 0..3 {
                    for j in 0..3 {
                        let dot: f64 = (0..3).map(|k| m[i][k] * m[j][k]).sum();
                        if (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1e-10 {
                            return Err(SpiceError::FormatParse(
                                "FK matrix is not orthonormal".into(),
                            ));
                        }
                    }
                }
                let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                if (det - 1.0).abs() > 1e-10 {
                    return Err(SpiceError::FormatParse(
                        "FK matrix must be right-handed".into(),
                    ));
                }
            }
            frames.push(FrameSpec {
                id,
                name,
                class,
                class_id,
                center_id,
                tk_spec,
            });
        }

        Ok(Self { frames })
    }

    /// Look up a frame by NAIF integer ID.
    pub fn frame_by_id(&self, id: i32) -> Option<&FrameSpec> {
        self.frames.iter().find(|frame| frame.id == id)
    }

    /// Look up a frame by name (case-insensitive).
    pub fn frame_by_name(&self, name: &str) -> Option<&FrameSpec> {
        let upper = name.to_ascii_uppercase();
        self.frames.iter().find(|frame| frame.name == upper)
    }
}

fn parse_frame_class(value: i32) -> Result<FrameClass, SpiceError> {
    match value {
        1 => Ok(FrameClass::Inertial),
        2 => Ok(FrameClass::Pck),
        3 => Ok(FrameClass::Ck),
        4 => Ok(FrameClass::TwoVector),
        5 => Ok(FrameClass::Dynamic),
        other => Err(SpiceError::FormatParse(format!(
            "unknown frame class {other}"
        ))),
    }
}

fn int_value(value: &TextValue, key: &str) -> Result<i32, SpiceError> {
    match value {
        TextValue::Integer(number) => i32::try_from(*number).map_err(|_| {
            SpiceError::FormatParse(format!("{key} integer {number} does not fit in i32"))
        }),
        TextValue::Float(number)
            if number.is_finite()
                && number.fract() == 0.0
                && *number >= i32::MIN as f64
                && *number <= i32::MAX as f64 =>
        {
            Ok(*number as i32)
        }
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be numeric, got {other:?}"
        ))),
    }
}

fn text_value<'a>(value: &'a TextValue, key: &str) -> Result<&'a str, SpiceError> {
    match value {
        TextValue::Text(text) => Ok(text.as_str()),
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be text, got {other:?}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::{FrameClass, FrameKernel, TkSpec};

    #[test]
    fn parse_minimal_inertial_frame() {
        let src = "\\begindata\nFRAME_1234_NAME = 'TEST'\nFRAME_1234_CLASS = 1\nFRAME_1234_CLASS_ID = 1234\nFRAME_1234_CENTER = 0\n";
        let kernel = FrameKernel::from_text(src).unwrap();
        let frame = kernel.frame_by_id(1234).unwrap();
        assert_eq!(frame.name, "TEST");
        assert_eq!(frame.class, FrameClass::Inertial);
    }

    #[test]
    fn parse_tk_matrix_frame() {
        let src = "\\begindata\nFRAME_2000_NAME = 'TEST_TK'\nFRAME_2000_CLASS = 4\nFRAME_2000_CLASS_ID = 2000\nFRAME_2000_CENTER = 399\nTKFRAME_2000_SPEC = 'MATRIX'\nTKFRAME_2000_RELATIVE = 'J2000'\nTKFRAME_2000_MATRIX = ( 1 0 0 0 1 0 0 0 1 )\n";
        let kernel = FrameKernel::from_text(src).unwrap();
        let frame = kernel.frame_by_name("test_tk").unwrap();
        match frame.tk_spec.as_ref().unwrap().0 {
            TkSpec::Matrix(matrix) => assert_eq!(matrix[0][0], 1.0),
        }
    }

    #[test]
    fn frame_lookup_by_id_and_name() {
        let src = "\\begindata\nFRAME_1234_NAME = 'TEST'\nFRAME_1234_CLASS = 1\nFRAME_1234_CLASS_ID = 1234\nFRAME_1234_CENTER = 0\n";
        let kernel = FrameKernel::from_text(src).unwrap();
        assert_eq!(kernel.frame_by_id(1234).unwrap().name, "TEST");
        assert_eq!(kernel.frame_by_name("test").unwrap().id, 1234);
    }
}
