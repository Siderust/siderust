// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Binary CK attitude-kernel support.
//!
//! ## Scientific scope
//!
//! This module provides discrete quaternion attitude lookup for spacecraft and
//! instrument frames stored in CK kernels.
//!
//! ## Technical scope
//!
//! V1 supports DAF/CK Type 1 segments only. Records are loaded into memory and
//! interpolated with SLERP.
//!
//! ## References
//!
//! - NAIF. *CK Required Reading*.
//! - NAIF. *DAF Required Reading*.

use std::path::Path;

use super::daf::DafRaw;
use super::SpiceError;

/// A single CK attitude record (Type 1).
#[derive(Debug, Clone)]
pub struct CkRecord {
    /// SCLK ticks (encoded as `f64` for comparison).
    pub sclkdp: f64,
    /// Unit quaternion `[q0, q1, q2, q3]` (scalar first).
    pub quaternion: [f64; 4],
}

/// A loaded CK Type 1 segment.
#[derive(Debug, Clone)]
pub struct CkSegment1 {
    /// Instrument NAIF ID (negative for spacecraft instruments).
    pub inst_id: i32,
    /// Reference frame ID.
    pub frame_id: i32,
    /// Coverage start (SCLK ticks encoded as `f64`).
    pub start_sclkdp: f64,
    /// Coverage end (SCLK ticks encoded as `f64`).
    pub end_sclkdp: f64,
    /// Records sorted by SCLK value.
    pub records: Vec<CkRecord>,
}

impl CkSegment1 {
    /// Interpolate attitude at SCLK ticks `sclkdp`.
    ///
    /// Uses SLERP between the two nearest bracketing records.
    pub fn evaluate(&self, sclkdp: f64) -> Result<[f64; 4], SpiceError> {
        let first = self.records.first().ok_or_else(|| SpiceError::Corrupted {
            message: format!("CK segment for instrument {} has no records", self.inst_id),
        })?;
        let last = self.records.last().expect("non-empty records");
        if sclkdp < first.sclkdp || sclkdp > last.sclkdp {
            return Err(SpiceError::OutOfCoverage {
                target: self.inst_id,
                center: self.frame_id,
                epoch_tdb_seconds: sclkdp,
                start_tdb_seconds: self.start_sclkdp,
                end_tdb_seconds: self.end_sclkdp,
            });
        }

        match self
            .records
            .binary_search_by(|record| record.sclkdp.total_cmp(&sclkdp))
        {
            Ok(index) => Ok(normalize(self.records[index].quaternion)),
            Err(insertion) => {
                if insertion == 0 || insertion >= self.records.len() {
                    return Err(SpiceError::OutOfCoverage {
                        target: self.inst_id,
                        center: self.frame_id,
                        epoch_tdb_seconds: sclkdp,
                        start_tdb_seconds: self.start_sclkdp,
                        end_tdb_seconds: self.end_sclkdp,
                    });
                }
                let left = &self.records[insertion - 1];
                let right = &self.records[insertion];
                let span = right.sclkdp - left.sclkdp;
                if span.abs() <= f64::EPSILON {
                    return Ok(normalize(left.quaternion));
                }
                let t = (sclkdp - left.sclkdp) / span;
                Ok(slerp(left.quaternion, right.quaternion, t))
            }
        }
    }
}

/// A loaded CK kernel ready to answer attitude queries.
#[derive(Debug)]
pub struct CkKernel {
    segments: Vec<CkSegment1>,
}

impl CkKernel {
    /// Parse a DAF/CK binary kernel from raw bytes.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self, SpiceError> {
        let daf = DafRaw::parse(bytes)?;
        if !daf
            .locator
            .trim()
            .to_ascii_uppercase()
            .starts_with("DAF/CK")
        {
            return Err(SpiceError::FormatParse(format!(
                "expected DAF/CK locator, got '{}'",
                daf.locator.trim()
            )));
        }

        let mut segments = Vec::new();
        for summary in &daf.raw_summaries {
            if summary.doubles.len() < 2 || summary.integers.len() < 6 {
                return Err(SpiceError::FormatParse(
                    "CK summary must contain ND=2 doubles and NI=6 integers".to_string(),
                ));
            }
            let inst_id = summary.integers[0];
            let frame_id = summary.integers[1];
            let data_type = summary.integers[2];
            let _av_flag = summary.integers[3];
            if data_type != 1 {
                continue;
            }
            segments.push(read_type1_segment(
                &daf,
                inst_id,
                frame_id,
                summary.doubles[0],
                summary.doubles[1],
                summary.start_word(),
                summary.end_word(),
            )?);
        }

        Ok(Self { segments })
    }

    /// Open a DAF/CK binary kernel from a filesystem path.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SpiceError> {
        let bytes = std::fs::read(path)?;
        Self::from_bytes(bytes)
    }

    /// Return the rotation quaternion for instrument `inst_id` at SCLK ticks `sclkdp`.
    pub fn rotation(&self, inst_id: i32, sclkdp: f64) -> Result<[f64; 4], SpiceError> {
        let mut coverage: Option<(f64, f64)> = None;
        for segment in &self.segments {
            if segment.inst_id != inst_id {
                continue;
            }
            coverage = Some(match coverage {
                Some((start, end)) => {
                    (start.min(segment.start_sclkdp), end.max(segment.end_sclkdp))
                }
                None => (segment.start_sclkdp, segment.end_sclkdp),
            });
            if sclkdp >= segment.start_sclkdp && sclkdp <= segment.end_sclkdp {
                return segment.evaluate(sclkdp);
            }
        }

        if let Some((start, end)) = coverage {
            return Err(SpiceError::OutOfCoverage {
                target: inst_id,
                center: 0,
                epoch_tdb_seconds: sclkdp,
                start_tdb_seconds: start,
                end_tdb_seconds: end,
            });
        }

        Err(SpiceError::UnsupportedKernelQuery {
            message: format!("no CK segment found for instrument {inst_id}"),
        })
    }

    /// All loaded segments.
    pub fn segments(&self) -> &[CkSegment1] {
        &self.segments
    }
}

fn read_type1_segment(
    daf: &DafRaw,
    inst_id: i32,
    frame_id: i32,
    start_sclkdp: f64,
    end_sclkdp: f64,
    start_word: usize,
    end_word: usize,
) -> Result<CkSegment1, SpiceError> {
    if end_word <= start_word {
        return Err(SpiceError::FormatParse(
            "CK segment has an invalid word range".to_string(),
        ));
    }
    let n_records_raw = daf.read_f64_at_word(end_word);
    if !n_records_raw.is_finite() || n_records_raw < 1.0 {
        return Err(SpiceError::FormatParse(format!(
            "invalid CK record count {n_records_raw}"
        )));
    }
    let n_records = n_records_raw as usize;
    let expected_words = n_records
        .checked_mul(5)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| SpiceError::FormatParse("CK record count overflow".to_string()))?;
    if end_word - start_word + 1 < expected_words {
        return Err(SpiceError::FormatParse(
            "CK segment word range is too small for its record count".to_string(),
        ));
    }

    let mut records = Vec::with_capacity(n_records);
    for index in 0..n_records {
        let word = start_word + index * 5;
        records.push(CkRecord {
            sclkdp: daf.read_f64_at_word(word),
            quaternion: normalize([
                daf.read_f64_at_word(word + 1),
                daf.read_f64_at_word(word + 2),
                daf.read_f64_at_word(word + 3),
                daf.read_f64_at_word(word + 4),
            ]),
        });
    }
    records.sort_by(|left, right| left.sclkdp.total_cmp(&right.sclkdp));

    Ok(CkSegment1 {
        inst_id,
        frame_id,
        start_sclkdp,
        end_sclkdp,
        records,
    })
}

fn normalize(quaternion: [f64; 4]) -> [f64; 4] {
    let norm = quaternion
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if norm <= f64::EPSILON {
        return quaternion;
    }
    [
        quaternion[0] / norm,
        quaternion[1] / norm,
        quaternion[2] / norm,
        quaternion[3] / norm,
    ]
}

fn slerp(left: [f64; 4], mut right: [f64; 4], t: f64) -> [f64; 4] {
    let left = normalize(left);
    right = normalize(right);
    let mut dot = left
        .iter()
        .zip(right.iter())
        .map(|(l, r)| l * r)
        .sum::<f64>();
    if dot < 0.0 {
        dot = -dot;
        right = [-right[0], -right[1], -right[2], -right[3]];
    }
    if dot >= 1.0 {
        return left;
    }
    if dot > 0.9995 {
        let interpolated = [
            left[0] + t * (right[0] - left[0]),
            left[1] + t * (right[1] - left[1]),
            left[2] + t * (right[2] - left[2]),
            left[3] + t * (right[3] - left[3]),
        ];
        return normalize(interpolated);
    }

    let theta = dot.acos();
    let sin_theta = theta.sin();
    let w_left = ((1.0 - t) * theta).sin() / sin_theta;
    let w_right = (t * theta).sin() / sin_theta;
    normalize([
        w_left * left[0] + w_right * right[0],
        w_left * left[1] + w_right * right[1],
        w_left * left[2] + w_right * right[2],
        w_left * left[3] + w_right * right[3],
    ])
}

#[cfg(test)]
mod tests {
    use super::CkKernel;

    fn write_i32_le(buf: &mut [u8], offset: usize, value: i32) {
        buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    }

    fn write_f64_le(buf: &mut [u8], offset: usize, value: f64) {
        buf[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    }

    fn synthetic_ck() -> Vec<u8> {
        let mut buf = vec![0u8; 3 * 1024];
        buf[0..8].copy_from_slice(b"DAF/CK  ");
        write_i32_le(&mut buf, 8, 2);
        write_i32_le(&mut buf, 12, 6);
        write_i32_le(&mut buf, 76, 2);

        let rec = &mut buf[1024..2048];
        write_f64_le(rec, 0, 0.0);
        write_f64_le(rec, 8, 0.0);
        write_f64_le(rec, 16, 1.0);

        let start_word = 257i32;
        let end_word = 267i32;
        write_f64_le(rec, 24, 0.0);
        write_f64_le(rec, 32, 10.0);
        write_i32_le(rec, 40, -82_000);
        write_i32_le(rec, 44, 1);
        write_i32_le(rec, 48, 1);
        write_i32_le(rec, 52, 0);
        write_i32_le(rec, 56, start_word);
        write_i32_le(rec, 60, end_word);

        let write_word = |word: i32, value: f64, bytes: &mut [u8]| {
            let offset = (word as usize - 1) * 8;
            bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
        };
        write_word(start_word, 0.0, &mut buf);
        write_word(start_word + 1, 1.0, &mut buf);
        write_word(start_word + 2, 0.0, &mut buf);
        write_word(start_word + 3, 0.0, &mut buf);
        write_word(start_word + 4, 0.0, &mut buf);
        write_word(start_word + 5, 10.0, &mut buf);
        write_word(start_word + 6, 0.0, &mut buf);
        write_word(start_word + 7, 1.0, &mut buf);
        write_word(start_word + 8, 0.0, &mut buf);
        write_word(start_word + 9, 0.0, &mut buf);
        write_word(end_word, 2.0, &mut buf);
        buf
    }

    #[test]
    fn parse_synthetic_ck_kernel() {
        let kernel = CkKernel::from_bytes(synthetic_ck()).unwrap();
        assert_eq!(kernel.segments().len(), 1);
        assert_eq!(kernel.segments()[0].inst_id, -82_000);
    }

    #[test]
    fn interpolate_type1_segment_with_slerp() {
        let kernel = CkKernel::from_bytes(synthetic_ck()).unwrap();
        let quaternion = kernel.rotation(-82_000, 5.0).unwrap();
        let half = std::f64::consts::FRAC_1_SQRT_2;
        assert!((quaternion[0] - half).abs() < 1.0e-12);
        assert!((quaternion[1] - half).abs() < 1.0e-12);
    }
}
