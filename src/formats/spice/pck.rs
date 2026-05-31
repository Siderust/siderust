// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Text-PCK body-orientation support.
//!
//! ## Scientific scope
//!
//! This module evaluates first-order IAU body pole and prime-meridian models
//! from text PCK kernels.
//!
//! ## Technical scope
//!
//! Only text PCK kernels are supported. V1 uses the leading linear terms of the
//! RA, Dec, and PM polynomials and ignores nutation terms.
//!
//! ## References
//!
//! - NAIF. *PCK Required Reading*.
//! - NAIF. *Kernel Required Reading*.

use std::collections::HashMap;

use super::text::{TextKernel, TextValue};
use super::SpiceError;

/// IAU orientation constants for one body.
#[derive(Debug, Clone)]
pub struct BodyOrientation {
    /// RA polynomial coefficients [a0, a1, ...] (degrees + degrees/century).
    pub pole_ra: Vec<f64>,
    /// Dec polynomial coefficients [d0, d1, ...] (degrees + degrees/century).
    pub pole_dec: Vec<f64>,
    /// PM polynomial coefficients [W0, W1, ...] (degrees + degrees/day).
    pub pm: Vec<f64>,
    /// Body radii in km: [equatorial, equatorial, polar].
    pub radii: Option<[f64; 3]>,
}

impl BodyOrientation {
    /// Rotation matrix from body-fixed frame to J2000 at epoch `tdb_s`.
    ///
    /// Uses first-order IAU pole RA/Dec/PM, ignoring nutation terms.
    pub fn rotation_to_j2000(&self, tdb_s: f64) -> [[f64; 3]; 3] {
        let t = tdb_s / (36_525.0 * 86_400.0);
        let d = tdb_s / 86_400.0;
        let ra_deg = self.pole_ra[0] + self.pole_ra.get(1).copied().unwrap_or(0.0) * t;
        let dec_deg = self.pole_dec[0] + self.pole_dec.get(1).copied().unwrap_or(0.0) * t;
        let w_deg = self.pm[0] + self.pm.get(1).copied().unwrap_or(0.0) * d;

        let ra = ra_deg.to_radians();
        let dec = dec_deg.to_radians();
        let w = w_deg.to_radians();

        mat3_mul(
            mat3_mul(
                r3(std::f64::consts::FRAC_PI_2 + ra),
                r1(std::f64::consts::FRAC_PI_2 - dec),
            ),
            r3(-w),
        )
    }
}

/// Parsed PCK text kernel containing body-orientation constants.
#[derive(Debug, Default, Clone)]
pub struct PckKernel {
    /// Parsed bodies keyed by NAIF ID.
    pub bodies: HashMap<i32, BodyOrientation>,
}

impl PckKernel {
    /// Parse a PCK text kernel.
    pub fn from_text(src: &str) -> Result<Self, SpiceError> {
        let kernel = TextKernel::parse(src)?;
        let mut bodies: HashMap<i32, BodyOrientation> = HashMap::new();

        for (key, value) in &kernel.data {
            let Some(rest) = key.strip_prefix("BODY") else {
                continue;
            };
            let Some((id_text, suffix)) = rest.split_once('_') else {
                continue;
            };
            let Ok(id) = id_text.parse::<i32>() else {
                continue;
            };
            let entry = bodies.entry(id).or_insert_with(|| BodyOrientation {
                pole_ra: Vec::new(),
                pole_dec: Vec::new(),
                pm: Vec::new(),
                radii: None,
            });
            match suffix {
                "POLE_RA" => entry.pole_ra = value_as_f64_array(value, key)?,
                "POLE_DEC" => entry.pole_dec = value_as_f64_array(value, key)?,
                "PM" => entry.pm = value_as_f64_array(value, key)?,
                "RADII" => entry.radii = Some(parse_radii(value, key)?),
                _ => {}
            }
        }

        Ok(Self { bodies })
    }

    /// Look up orientation constants for a NAIF body ID.
    pub fn body(&self, id: i32) -> Option<&BodyOrientation> {
        self.bodies.get(&id)
    }
}

fn value_as_f64_array(value: &TextValue, key: &str) -> Result<Vec<f64>, SpiceError> {
    match value {
        TextValue::Integer(number) => Ok(vec![*number as f64]),
        TextValue::Float(number) => Ok(vec![*number]),
        TextValue::Array(values) => values
            .iter()
            .map(|entry| match entry {
                TextValue::Integer(number) => Ok(*number as f64),
                TextValue::Float(number) => Ok(*number),
                other => Err(SpiceError::FormatParse(format!(
                    "{key} must contain only numeric entries, got {other:?}"
                ))),
            })
            .collect(),
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be numeric, got {other:?}"
        ))),
    }
}

fn parse_radii(value: &TextValue, key: &str) -> Result<[f64; 3], SpiceError> {
    let values = value_as_f64_array(value, key)?;
    match values.as_slice() {
        [equatorial, polar] => Ok([*equatorial, *equatorial, *polar]),
        [x, y, z] => Ok([*x, *y, *z]),
        _ => Err(SpiceError::FormatParse(format!(
            "{key} must contain 2 or 3 radii values"
        ))),
    }
}

fn r1(angle: f64) -> [[f64; 3]; 3] {
    let (sin, cos) = angle.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, cos, sin], [0.0, -sin, cos]]
}

fn r3(angle: f64) -> [[f64; 3]; 3] {
    let (sin, cos) = angle.sin_cos();
    [[cos, sin, 0.0], [-sin, cos, 0.0], [0.0, 0.0, 1.0]]
}

fn mat3_mul(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += left[i][k] * right[k][j];
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::PckKernel;

    #[test]
    fn parse_minimal_pck_with_earth_data() {
        let src = "\\begindata\nBODY399_POLE_RA = ( 0.0 0.0 )\nBODY399_POLE_DEC = ( 90.0 0.0 )\nBODY399_PM = ( 190.0 360.9856235 )\nBODY399_RADII = ( 6378.137 6378.137 6356.752 )\n";
        let kernel = PckKernel::from_text(src).unwrap();
        assert!(kernel.body(399).is_some());
    }

    #[test]
    fn rotation_is_orthonormal_at_j2000() {
        let src = "\\begindata\nBODY399_POLE_RA = ( 0.0 0.0 )\nBODY399_POLE_DEC = ( 90.0 0.0 )\nBODY399_PM = ( 0.0 1.0 )\n";
        let kernel = PckKernel::from_text(src).unwrap();
        let matrix = kernel.body(399).unwrap().rotation_to_j2000(0.0);
        let det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        let row0_norm = (matrix[0][0].powi(2) + matrix[0][1].powi(2) + matrix[0][2].powi(2)).sqrt();
        assert!((row0_norm - 1.0).abs() < 1.0e-12);
        assert!((det - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn body_lookup_returns_correct_entry() {
        let src = "\\begindata\nBODY301_POLE_RA = ( 269.9949 0.0031 )\nBODY301_POLE_DEC = ( 66.5392 0.0130 )\nBODY301_PM = ( 38.3213 13.17635815 )\n";
        let kernel = PckKernel::from_text(src).unwrap();
        assert_eq!(kernel.body(301).unwrap().pole_ra[0], 269.9949);
    }
}
