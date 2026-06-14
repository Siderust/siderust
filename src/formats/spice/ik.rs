// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Instrument-kernel support.
//!
//! ## Scientific scope
//!
//! This module parses the minimal instrument field-of-view metadata needed for
//! SPICE V1 geometry context assembly.
//!
//! ## Technical scope
//!
//! Only text IK kernels are supported. V1 extracts FOV shape, frame, boresight,
//! and angle fields for each instrument ID encountered.
//!
//! ## References
//!
//! - NAIF. *IK Required Reading*.
//! - NAIF. *Kernel Required Reading*.

use super::text::{TextKernel, TextValue};
use super::SpiceError;

/// Instrument field-of-view shape.
#[derive(Debug, Clone)]
pub enum FovShape {
    /// Circular field of view.
    Circle,
    /// Rectangular field of view.
    Rectangle,
    /// Elliptical field of view.
    Ellipse,
    /// Polygonal field of view.
    Polygon,
    /// Any other shape string.
    Other(String),
}

/// Instrument geometric parameters parsed from an IK text kernel.
#[derive(Debug, Clone)]
pub struct InstrumentKernel {
    /// NAIF instrument ID (negative integer).
    pub inst_id: i32,
    /// FOV shape.
    pub fov_shape: Option<FovShape>,
    /// Frame name for the FOV.
    pub fov_frame: Option<String>,
    /// Boresight vector (body-frame).
    pub boresight: Option<[f64; 3]>,
    /// Half-angle in degrees (for circle/rectangle).
    pub fov_ref_angle: Option<f64>,
    /// Cross-angle in degrees (for rectangle only).
    pub fov_cross_angle: Option<f64>,
}

/// Collection of instrument specs from one IK kernel.
#[derive(Debug, Default, Clone)]
pub struct IkKernel {
    /// Parsed instrument definitions.
    pub instruments: Vec<InstrumentKernel>,
}

impl IkKernel {
    /// Parse an IK text kernel.
    pub fn from_text(src: &str) -> Result<Self, SpiceError> {
        let kernel = TextKernel::parse(src)?;
        let mut ids = Vec::new();
        for key in kernel.data.keys() {
            let Some(rest) = key.strip_prefix("INS") else {
                continue;
            };
            let Some((id_text, _suffix)) = rest.split_once('_') else {
                continue;
            };
            let Ok(id) = id_text.parse::<i32>() else {
                continue;
            };
            ids.push(id);
        }
        ids.sort_unstable();
        ids.dedup();

        let mut instruments = Vec::with_capacity(ids.len());
        for inst_id in ids {
            let shape = kernel
                .get(&format!("INS{inst_id}_FOV_SHAPE"))
                .map(|value| parse_shape(value, inst_id))
                .transpose()?;
            let frame = kernel
                .get(&format!("INS{inst_id}_FOV_FRAME"))
                .map(|value| {
                    text_value(value, &format!("INS{inst_id}_FOV_FRAME"))
                        .map(|text| text.to_ascii_uppercase())
                })
                .transpose()?;
            let boresight = kernel
                .get_f64_array(&format!("INS{inst_id}_BORESIGHT"))
                .map(|values| {
                    if values.len() == 3 {
                        Ok([values[0], values[1], values[2]])
                    } else {
                        Err(SpiceError::FormatParse(format!(
                            "INS{inst_id}_BORESIGHT must contain 3 values"
                        )))
                    }
                })
                .transpose()?;
            let fov_ref_angle = kernel
                .get(&format!("INS{inst_id}_FOV_REF_ANGLE"))
                .map(|value| scalar_value(value, &format!("INS{inst_id}_FOV_REF_ANGLE")))
                .transpose()?;
            let fov_cross_angle = kernel
                .get(&format!("INS{inst_id}_FOV_CROSS_ANGLE"))
                .map(|value| scalar_value(value, &format!("INS{inst_id}_FOV_CROSS_ANGLE")))
                .transpose()?;

            instruments.push(InstrumentKernel {
                inst_id,
                fov_shape: shape,
                fov_frame: frame,
                boresight,
                fov_ref_angle,
                fov_cross_angle,
            });
        }

        Ok(Self { instruments })
    }

    /// Look up an instrument by NAIF ID.
    pub fn instrument(&self, id: i32) -> Option<&InstrumentKernel> {
        self.instruments
            .iter()
            .find(|instrument| instrument.inst_id == id)
    }
}

fn parse_shape(value: &TextValue, inst_id: i32) -> Result<FovShape, SpiceError> {
    let text = text_value(value, &format!("INS{inst_id}_FOV_SHAPE"))?.to_ascii_uppercase();
    Ok(match text.as_str() {
        "CIRCLE" => FovShape::Circle,
        "RECTANGLE" => FovShape::Rectangle,
        "ELLIPSE" => FovShape::Ellipse,
        "POLYGON" => FovShape::Polygon,
        other => FovShape::Other(other.to_string()),
    })
}

fn text_value<'a>(value: &'a TextValue, key: &str) -> Result<&'a str, SpiceError> {
    match value {
        TextValue::Text(text) => Ok(text.as_str()),
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be text, got {other:?}"
        ))),
    }
}

fn scalar_value(value: &TextValue, key: &str) -> Result<f64, SpiceError> {
    match value {
        TextValue::Integer(number) => Ok(*number as f64),
        TextValue::Float(number) => Ok(*number),
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be numeric, got {other:?}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::{FovShape, IkKernel};

    #[test]
    fn parse_minimal_ik_kernel() {
        let src = "\\begindata\nINS-82000_FOV_SHAPE = 'RECTANGLE'\nINS-82000_FOV_FRAME = 'CASSINI_ISS_NAC'\nINS-82000_BORESIGHT = ( 0 0 1 )\nINS-82000_FOV_REF_ANGLE = 0.35\nINS-82000_FOV_CROSS_ANGLE = 0.35\n";
        let kernel = IkKernel::from_text(src).unwrap();
        let instrument = kernel.instrument(-82_000).unwrap();
        assert!(matches!(instrument.fov_shape, Some(FovShape::Rectangle)));
        assert_eq!(instrument.fov_frame.as_deref(), Some("CASSINI_ISS_NAC"));
        assert_eq!(instrument.boresight, Some([0.0, 0.0, 1.0]));
    }
}
