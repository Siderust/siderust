// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spacecraft-clock text-kernel support.
//!
//! ## Scientific scope
//!
//! This module provides the minimal Type-1 SCLK mapping needed to convert
//! between SPICE ephemeris time (TDB seconds past J2000) and encoded spacecraft
//! clock ticks.
//!
//! ## Technical scope
//!
//! Only text SCLK kernels with `SCLK01_COEFFICIENTS_<id>` entries are supported.
//! Conversion is piecewise linear within each partition.
//!
//! ## References
//!
//! - NAIF. *SCLK Required Reading*.
//! - NAIF. *Kernel Required Reading*.

use super::text::{TextKernel, TextValue};
use super::SpiceError;

/// One coefficient entry for an SCLK Type 1 clock.
#[derive(Debug, Clone, Copy)]
pub struct SclkCoeff {
    /// Encoded SCLK ticks at partition start.
    pub sclkdp: f64,
    /// TDB seconds past J2000 at partition start.
    pub et: f64,
    /// SCLK ticks per TDB second.
    pub rate: f64,
}

/// Parsed SCLK kernel for one spacecraft clock.
#[derive(Debug, Clone)]
pub struct SclkKernel {
    /// Spacecraft NAIF ID (negative integer, e.g. -82 for Cassini).
    pub sc_id: i32,
    /// Coefficient table, sorted by `et` ascending.
    pub coefficients: Vec<SclkCoeff>,
}

impl SclkKernel {
    /// Parse a SCLK text kernel. Returns one kernel per spacecraft ID found.
    pub fn from_text(src: &str) -> Result<Vec<Self>, SpiceError> {
        let kernel = TextKernel::parse(src)?;
        let mut ids = Vec::new();
        for key in kernel.data.keys() {
            if let Some(id) = key
                .strip_prefix("SCLK01_COEFFICIENTS_")
                .and_then(|rest| rest.parse::<i32>().ok())
            {
                ids.push(id);
            }
        }
        ids.sort_unstable();
        ids.dedup();

        let mut out = Vec::with_capacity(ids.len());
        for sc_id in ids {
            let key = format!("SCLK01_COEFFICIENTS_{sc_id}");
            let values = value_as_f64_array(
                kernel
                    .get(&key)
                    .ok_or_else(|| SpiceError::FormatParse(format!("missing {key}")))?,
                &key,
            )?;
            if values.len() % 3 != 0 {
                return Err(SpiceError::FormatParse(format!(
                    "{key} must contain triples of (sclkdp, et, rate)"
                )));
            }
            let mut coefficients = Vec::with_capacity(values.len() / 3);
            for triple in values.chunks(3) {
                coefficients.push(SclkCoeff {
                    sclkdp: triple[0],
                    et: triple[1],
                    rate: triple[2],
                });
            }
            coefficients.sort_by(|left, right| left.et.total_cmp(&right.et));
            out.push(Self {
                sc_id,
                coefficients,
            });
        }
        Ok(out)
    }

    /// Convert TDB seconds past J2000 to encoded SCLK ticks.
    pub fn et_to_sclkdp(&self, et: f64) -> Result<f64, SpiceError> {
        let coeff = self.partition_for_et(et)?;
        Ok(coeff.sclkdp + (et - coeff.et) * coeff.rate)
    }

    /// Convert encoded SCLK ticks to TDB seconds past J2000.
    pub fn sclkdp_to_et(&self, sclkdp: f64) -> Result<f64, SpiceError> {
        let coeff = self.partition_for_sclkdp(sclkdp)?;
        Ok(coeff.et + (sclkdp - coeff.sclkdp) / coeff.rate)
    }

    fn partition_for_et(&self, et: f64) -> Result<SclkCoeff, SpiceError> {
        let first =
            self.coefficients
                .first()
                .copied()
                .ok_or_else(|| SpiceError::TimeConversion {
                    message: format!("SCLK {} has no coefficients", self.sc_id),
                })?;
        if et < first.et {
            return Err(SpiceError::TimeConversion {
                message: format!("epoch {et} is before SCLK {} coverage", self.sc_id),
            });
        }
        for window in self.coefficients.windows(2) {
            let current = window[0];
            let next = window[1];
            if et >= current.et && et < next.et {
                return Ok(current);
            }
        }
        Ok(*self.coefficients.last().expect("non-empty coefficients"))
    }

    fn partition_for_sclkdp(&self, sclkdp: f64) -> Result<SclkCoeff, SpiceError> {
        let first =
            self.coefficients
                .first()
                .copied()
                .ok_or_else(|| SpiceError::TimeConversion {
                    message: format!("SCLK {} has no coefficients", self.sc_id),
                })?;
        if sclkdp < first.sclkdp {
            return Err(SpiceError::TimeConversion {
                message: format!(
                    "encoded ticks {sclkdp} are before SCLK {} coverage",
                    self.sc_id
                ),
            });
        }
        for window in self.coefficients.windows(2) {
            let current = window[0];
            let next = window[1];
            if sclkdp >= current.sclkdp && sclkdp < next.sclkdp {
                return Ok(current);
            }
        }
        Ok(*self.coefficients.last().expect("non-empty coefficients"))
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

#[cfg(test)]
mod tests {
    use super::SclkKernel;

    #[test]
    fn parse_minimal_sclk_with_two_partitions() {
        let src = "\\begindata\nSCLK01_COEFFICIENTS_-82 = ( 0 0 2 100 50 4 )\n";
        let kernels = SclkKernel::from_text(src).unwrap();
        assert_eq!(kernels.len(), 1);
        assert_eq!(kernels[0].sc_id, -82);
        assert_eq!(kernels[0].coefficients.len(), 2);
    }

    #[test]
    fn et_and_sclkdp_are_approximate_inverses() {
        let src = "\\begindata\nSCLK01_COEFFICIENTS_-82 = ( 0 0 2 100 50 4 )\n";
        let mut kernels = SclkKernel::from_text(src).unwrap();
        let kernel = kernels.remove(0);
        let et = 25.0;
        let sclkdp = kernel.et_to_sclkdp(et).unwrap();
        let round_trip = kernel.sclkdp_to_et(sclkdp).unwrap();
        assert!((round_trip - et).abs() < 1.0e-12);
    }
}
