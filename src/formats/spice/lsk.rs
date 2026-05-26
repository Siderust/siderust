// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Leapseconds-kernel support.
//!
//! ## Scientific scope
//!
//! This module provides the minimal leapseconds support needed to map between
//! UTC-like offsets and TDB-style J2000 seconds for SPICE V1 use cases.
//!
//! ## Technical scope
//!
//! Only text LSK kernels are supported. The parser reads `DELTET/DELTA_T_A`
//! and `DELTET/DELTA_AT`, stores a sorted leap-second table, and exposes
//! simple offset queries.
//!
//! ## References
//!
//! - NAIF. *Leapseconds Kernel Required Reading*.
//! - NAIF. *Kernel Required Reading*.

use super::text::{TextKernel, TextValue};
use super::SpiceError;

/// Parsed leapseconds kernel.
#[derive(Debug, Clone)]
pub struct LeapSecondKernel {
    /// Constant TDT−TAI offset in seconds (DELTET/DELTA_T_A, ≈ 32.184 s).
    pub delta_t_a: f64,
    /// Sorted leapseconds table: `(epoch_j2000_tdb_seconds, tai_minus_utc)`.
    pub leap_table: Vec<(f64, f64)>,
}

impl LeapSecondKernel {
    /// Parse an LSK text kernel.
    pub fn from_text(src: &str) -> Result<Self, SpiceError> {
        let kernel = TextKernel::parse(src)?;
        let delta_t_a = scalar_as_f64(
            kernel.get("DELTET/DELTA_T_A").ok_or_else(|| {
                SpiceError::FormatParse("LSK missing DELTET/DELTA_T_A".to_string())
            })?,
            "DELTET/DELTA_T_A",
        )?;
        let delta_at = match kernel
            .get("DELTET/DELTA_AT")
            .ok_or_else(|| SpiceError::FormatParse("LSK missing DELTET/DELTA_AT".to_string()))?
        {
            TextValue::Array(values) => values,
            _ => {
                return Err(SpiceError::FormatParse(
                    "DELTET/DELTA_AT must be an array".to_string(),
                ));
            }
        };

        if delta_at.len() % 2 != 0 {
            return Err(SpiceError::FormatParse(
                "DELTET/DELTA_AT must contain delta/epoch pairs".to_string(),
            ));
        }

        let mut leap_table = Vec::with_capacity(delta_at.len() / 2);
        for pair in delta_at.chunks(2) {
            let delta = scalar_as_f64(&pair[0], "DELTET/DELTA_AT delta")?;
            let epoch = epoch_as_j2000_seconds(&pair[1])?;
            leap_table.push((epoch, delta));
        }
        leap_table.sort_by(|left, right| left.0.total_cmp(&right.0));

        Ok(Self {
            delta_t_a,
            leap_table,
        })
    }

    /// Convert TDB seconds past J2000 to the TAI-UTC offset valid at that epoch.
    pub fn tai_minus_utc_at(&self, tdb_seconds: f64) -> f64 {
        let mut best = 0.0;
        for &(epoch, delta) in &self.leap_table {
            if epoch <= tdb_seconds {
                best = delta;
            } else {
                break;
            }
        }
        best
    }

    /// TDB seconds past J2000 to total seconds ahead of UTC at that epoch.
    pub fn tdb_minus_utc_at(&self, tdb_seconds: f64) -> f64 {
        self.delta_t_a + self.tai_minus_utc_at(tdb_seconds)
    }
}

fn scalar_as_f64(value: &TextValue, key: &str) -> Result<f64, SpiceError> {
    match value {
        TextValue::Integer(number) => Ok(*number as f64),
        TextValue::Float(number) => Ok(*number),
        other => Err(SpiceError::FormatParse(format!(
            "{key} must be numeric, got {other:?}"
        ))),
    }
}

fn epoch_as_j2000_seconds(value: &TextValue) -> Result<f64, SpiceError> {
    match value {
        TextValue::Integer(number) => Ok(*number as f64),
        TextValue::Float(number) => Ok(*number),
        TextValue::Text(text) => parse_epoch_text(text),
        TextValue::Array(_) => Err(SpiceError::FormatParse(
            "LSK epoch entries cannot be arrays".to_string(),
        )),
    }
}

fn parse_epoch_text(text: &str) -> Result<f64, SpiceError> {
    let trimmed = text.trim();
    if !trimmed.starts_with('@') {
        return Err(SpiceError::FormatParse(format!(
            "unsupported LSK epoch token '{trimmed}'"
        )));
    }
    let date = &trimmed[1..];
    let mut parts = date.split('-');
    let year = parts
        .next()
        .ok_or_else(|| SpiceError::FormatParse(format!("bad LSK epoch '{trimmed}'")))?
        .parse::<i32>()
        .map_err(|_| SpiceError::FormatParse(format!("bad LSK year in '{trimmed}'")))?;
    let month = month_number(
        parts
            .next()
            .ok_or_else(|| SpiceError::FormatParse(format!("bad LSK epoch '{trimmed}'")))?,
    )?;
    let day_text = parts
        .next()
        .ok_or_else(|| SpiceError::FormatParse(format!("bad LSK epoch '{trimmed}'")))?;
    let day = day_text
        .parse::<i32>()
        .map_err(|_| SpiceError::FormatParse(format!("bad LSK day in '{trimmed}'")))?;

    let jd = gregorian_to_jd(year, month, day);
    Ok((jd - 2_451_545.0) * 86_400.0)
}

fn month_number(month: &str) -> Result<i32, SpiceError> {
    match month.trim().to_ascii_uppercase().as_str() {
        "JAN" => Ok(1),
        "FEB" => Ok(2),
        "MAR" => Ok(3),
        "APR" => Ok(4),
        "MAY" => Ok(5),
        "JUN" => Ok(6),
        "JUL" => Ok(7),
        "AUG" => Ok(8),
        "SEP" => Ok(9),
        "OCT" => Ok(10),
        "NOV" => Ok(11),
        "DEC" => Ok(12),
        other => Err(SpiceError::FormatParse(format!(
            "unknown LSK month abbreviation '{other}'"
        ))),
    }
}

fn gregorian_to_jd(year: i32, month: i32, day: i32) -> f64 {
    let term1 = 367 * year;
    let term2 = 7 * (year + (month + 9) / 12) / 4;
    let term3 = 275 * month / 9;
    (term1 - term2 + term3 + day + 1_721_013) as f64 + 0.5
}

#[cfg(test)]
mod tests {
    use super::LeapSecondKernel;

    #[test]
    fn parse_minimal_lsk_with_one_entry() {
        let src = "\\begindata\nDELTET/DELTA_T_A = 32.184\nDELTET/DELTA_AT = ( 37 @2017-JAN-1 )\n";
        let kernel = LeapSecondKernel::from_text(src).unwrap();
        assert_eq!(kernel.delta_t_a, 32.184);
        assert_eq!(kernel.leap_table.len(), 1);
        assert_eq!(kernel.leap_table[0].1, 37.0);
    }

    #[test]
    fn tai_minus_utc_returns_offset_before_and_after_epoch() {
        let src = "\\begindata\nDELTET/DELTA_T_A = 32.184\nDELTET/DELTA_AT = ( 37 @2017-JAN-1 )\n";
        let kernel = LeapSecondKernel::from_text(src).unwrap();
        let epoch = kernel.leap_table[0].0;
        assert_eq!(kernel.tai_minus_utc_at(epoch - 1.0), 0.0);
        assert_eq!(kernel.tai_minus_utc_at(epoch + 1.0), 37.0);
    }

    #[test]
    fn tdb_minus_utc_adds_delta_t_a() {
        let src = "\\begindata\nDELTET/DELTA_T_A = 32.184\nDELTET/DELTA_AT = ( 37 @2017-JAN-1 )\n";
        let kernel = LeapSecondKernel::from_text(src).unwrap();
        let epoch = kernel.leap_table[0].0 + 1.0;
        assert!((kernel.tdb_minus_utc_at(epoch) - 69.184).abs() < 1.0e-12);
    }
}
