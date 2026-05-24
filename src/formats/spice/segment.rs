// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SPK segment decoding and Chebyshev evaluation.
//!
//! ## Technical scope
//!
//! This module decodes the SPK segment payloads used by the low-level kernel
//! reader. It stays at the NAIF wire-format level: epochs are raw TDB seconds
//! past J2000 and evaluations return raw `[f64; 6]` state vectors in km and
//! km/s.
//!
//! ## API contract
//!
//! [`SpkSegment::evaluate`] is intentionally a low-level wire-format entry
//! point. Higher-level scientific callers should prefer typed ephemeris
//! providers instead of working with raw kernel records directly.
//!
//! This module covers the two SPK data types relevant to the JPL DE
//! planetary kernels exercised by POD:
//!
//! * **Type 2** — Chebyshev polynomials for position only; velocity is
//!   recovered analytically from the position polynomial derivative.
//! * **Type 3** — Chebyshev polynomials for position **and** velocity,
//!   with independent coefficient sets.
//!
//! Higher SPK types (9 / 13 / 21 / etc.) are out of scope; loading is
//! still permitted, but evaluation returns
//! [`SpiceError::UnsupportedDataType`].
//!
//! # Per-record layout
//!
//! Both Type 2 and Type 3 segments are stored as a packed sequence of
//! fixed-size records. Each record contains the midpoint and radius
//! (TDB seconds) of the interval it covers, followed by the polynomial
//! coefficients:
//!
//! ```text
//! [ MID, RADIUS, c0_x, c1_x, ..., c0_y, c1_y, ..., c0_z, c1_z, ...]
//!                  ^-- Type 2: 3 components (x, y, z)
//! [ MID, RADIUS, c0_x, ..., c0_y, ..., c0_z, ..., c0_vx, ..., c0_vy, ..., c0_vz, ...]
//!                  ^-- Type 3: 6 components
//! ```
//!
//! The polynomial argument is `tau = (et - MID) / RADIUS ∈ [-1, 1]`.
//! Velocity for Type 2 is `dPos/dtau / RADIUS` (in km/s).

use super::error::SpiceError;
use crate::formats::spice::daf::{Daf, Summary};

/// A typed SPK segment. Currently restricted to Chebyshev variants.
#[derive(Debug, Clone)]
pub enum SpkSegment {
    /// SPK Type 2 — Chebyshev polynomial for position; velocity is
    /// derived analytically.
    Type2(ChebSegment),
    /// SPK Type 3 — Chebyshev polynomials for position and velocity
    /// stored side by side (no analytic differentiation needed).
    Type3(ChebSegment),
}

impl SpkSegment {
    /// Coverage start (TDB seconds past J2000) of the segment.
    pub fn start_tdb_seconds(&self) -> f64 {
        self.cheb().init
    }

    /// Coverage end (TDB seconds past J2000) of the segment.
    pub fn end_tdb_seconds(&self) -> f64 {
        let c = self.cheb();
        c.init + c.intlen * c.n_records as f64
    }

    /// Borrow the underlying Chebyshev segment data.
    pub fn cheb(&self) -> &ChebSegment {
        match self {
            SpkSegment::Type2(c) | SpkSegment::Type3(c) => c,
        }
    }

    /// Evaluate the segment at `et` (TDB seconds past J2000).
    ///
    /// Returns `[x, y, z, vx, vy, vz]` with positions in **kilometers**
    /// and velocities in **kilometers per second**, expressed in the
    /// segment's reference frame (typically J2000 ≡ ICRS for DE
    /// kernels).
    ///
    /// # Safety / Wire-format API
    ///
    /// This method intentionally exposes raw SPK record values. Prefer a typed
    /// [`crate::pod::providers::EphemerisProvider`] implementation when you
    /// need the scientific API.
    pub fn evaluate(&self, et: f64) -> Result<[f64; 6], SpiceError> {
        match self {
            SpkSegment::Type2(c) => c.evaluate_type2(et),
            SpkSegment::Type3(c) => c.evaluate_type3(et),
        }
    }
}

/// Chebyshev segment data shared by Type 2 and Type 3.
#[derive(Debug, Clone)]
pub struct ChebSegment {
    /// Initial epoch of the segment (TDB seconds past J2000).
    pub init: f64,
    /// Length (seconds) of each record's coverage interval.
    pub intlen: f64,
    /// Doubles per record (`2 + components * ncoeff`).
    pub rsize: usize,
    /// Number of coefficients per Chebyshev series.
    pub ncoeff: usize,
    /// Number of records in the segment.
    pub n_records: usize,
    /// Components per record (`3` for Type 2, `6` for Type 3).
    pub components: u8,
    /// Flat record data (`n_records * rsize` doubles).
    pub records: Vec<f64>,
}

impl ChebSegment {
    /// Locate the record covering `et` and return `(record_index, mid, radius)`.
    fn locate(&self, et: f64) -> Result<(usize, f64, f64), SpiceError> {
        if !et.is_finite() {
            return Err(SpiceError::Corrupted {
                message: format!("epoch is not finite: {et}"),
            });
        }
        let end = self.init + self.intlen * self.n_records as f64;
        if et < self.init || et > end {
            return Err(SpiceError::OutOfCoverage {
                target: 0,
                center: 0,
                epoch_tdb_seconds: et,
                start_tdb_seconds: self.init,
                end_tdb_seconds: end,
            });
        }
        let mut idx = ((et - self.init) / self.intlen).floor() as i64;
        if idx < 0 {
            idx = 0;
        }
        if idx as usize >= self.n_records {
            idx = self.n_records as i64 - 1;
        }
        let off = (idx as usize) * self.rsize;
        let mid = self.records[off];
        let radius = self.records[off + 1];
        if !radius.is_finite() || radius <= 0.0 {
            return Err(SpiceError::Corrupted {
                message: format!("record {idx} has non-positive radius {radius} (mid={mid})"),
            });
        }
        Ok((idx as usize, mid, radius))
    }

    /// Evaluate as a Type 2 (Chebyshev position) segment.
    pub fn evaluate_type2(&self, et: f64) -> Result<[f64; 6], SpiceError> {
        debug_assert_eq!(self.components, 3);
        let (idx, mid, radius) = self.locate(et)?;
        let off = idx * self.rsize + 2;
        let tau = (et - mid) / radius;
        let n = self.ncoeff;
        let mut out = [0.0_f64; 6];
        for i in 0..3 {
            let coeffs = &self.records[off + i * n..off + (i + 1) * n];
            let (p, dp) = cheby::core::evaluate_both(coeffs, tau);
            out[i] = p;
            out[i + 3] = dp / radius;
        }
        Ok(out)
    }

    /// Evaluate as a Type 3 (Chebyshev position+velocity) segment.
    pub fn evaluate_type3(&self, et: f64) -> Result<[f64; 6], SpiceError> {
        debug_assert_eq!(self.components, 6);
        let (idx, mid, radius) = self.locate(et)?;
        let off = idx * self.rsize + 2;
        let tau = (et - mid) / radius;
        let n = self.ncoeff;
        let mut out = [0.0_f64; 6];
        for (i, slot) in out.iter_mut().enumerate() {
            let coeffs = &self.records[off + i * n..off + (i + 1) * n];
            *slot = cheby::core::evaluate(coeffs, tau);
        }
        Ok(out)
    }
}

/// Read an SPK segment summary into a typed [`SpkSegment`].
///
/// Supports SPK Type 2 and Type 3; any other type is rejected with
/// [`SpiceError::UnsupportedDataType`] so the caller can index the
/// kernel without losing information about what is and isn't loadable.
///
/// # Examples
///
/// ```rust
/// use siderust::formats::spice::{daf::{Daf, Summary}, segment_for_summary, SpiceError};
/// // Build a synthetic single-record Type 2 segment with one
/// // coefficient per axis (a constant polynomial).
/// // rsize = 2 (mid,radius) + 3 components * 1 coeff = 5 doubles.
/// let mid = 0.5_f64;
/// let radius = 0.5_f64;
/// let coeffs_x = 100.0_f64;
/// let coeffs_y = 200.0_f64;
/// let coeffs_z = 300.0_f64;
/// let mut buf = vec![0u8; 9 * 8];
/// let mut write = |word: usize, v: f64| {
///     buf[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
/// };
/// // Records (words 1..=5)
/// write(1, mid);
/// write(2, radius);
/// write(3, coeffs_x);
/// write(4, coeffs_y);
/// write(5, coeffs_z);
/// // Trailer (init, intlen, rsize, n_records) at words 6..=9
/// write(6, 0.0);   // init
/// write(7, 1.0);   // intlen
/// write(8, 5.0);   // rsize
/// write(9, 1.0);   // n_records
///
/// let daf = Daf { nd: 2, ni: 6, summaries: vec![] };
/// let summary = Summary {
///     start_et: 0.0, end_et: 1.0,
///     target_id: 399, center_id: 0, frame_id: 1, data_type: 2,
///     start_word: 1, end_word: 9,
/// };
/// let segment = segment_for_summary(&buf, &daf, &summary)?;
/// let state = segment.evaluate(0.5)?;
/// assert!((state[0] - 100.0).abs() < 1e-12);
/// assert!((state[1] - 200.0).abs() < 1e-12);
/// assert!((state[2] - 300.0).abs() < 1e-12);
/// # Ok::<_, SpiceError>(())
/// ```
pub fn segment_for_summary(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
) -> Result<SpkSegment, SpiceError> {
    match summary.data_type {
        2 => Ok(SpkSegment::Type2(read_chebyshev(
            file_data, daf, summary, 3,
        )?)),
        3 => Ok(SpkSegment::Type3(read_chebyshev(
            file_data, daf, summary, 6,
        )?)),
        other => Err(SpiceError::UnsupportedDataType { data_type: other }),
    }
}

fn read_chebyshev(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
    components: u8,
) -> Result<ChebSegment, SpiceError> {
    let end = summary.end_word;
    if end < 4 {
        return Err(SpiceError::Parse {
            message: format!("segment trailer at end_word={end} is too small"),
        });
    }
    let n_records = daf.read_f64_at_word(file_data, end) as usize;
    let rsize = daf.read_f64_at_word(file_data, end - 1) as usize;
    let intlen = daf.read_f64_at_word(file_data, end - 2);
    let init = daf.read_f64_at_word(file_data, end - 3);

    if !(5..=4096).contains(&rsize) {
        return Err(SpiceError::Parse {
            message: format!("implausible rsize={rsize} for SPK Chebyshev segment"),
        });
    }
    let comp = components as usize;
    if !(rsize - 2).is_multiple_of(comp) {
        return Err(SpiceError::Parse {
            message: format!("rsize={rsize} is not 2 + {comp}*k for any k (components={comp})"),
        });
    }
    let ncoeff = (rsize - 2) / comp;
    if n_records == 0 || n_records > 10_000_000 {
        return Err(SpiceError::Parse {
            message: format!("implausible n_records={n_records}"),
        });
    }
    if !(intlen.is_finite() && intlen > 0.0) {
        return Err(SpiceError::Parse {
            message: format!("non-positive intlen={intlen}"),
        });
    }
    if !init.is_finite() {
        return Err(SpiceError::Parse {
            message: format!("non-finite init={init}"),
        });
    }

    let total = n_records
        .checked_mul(rsize)
        .ok_or_else(|| SpiceError::Parse {
            message: format!("segment record count overflow: n_records={n_records} rsize={rsize}"),
        })?;
    let data_start_word = summary.start_word;
    if data_start_word + total - 1 > summary.end_word {
        return Err(SpiceError::Parse {
            message: format!(
                "segment data window [{}..{}] does not fit in summary [{}..{}]",
                data_start_word,
                data_start_word + total - 1,
                summary.start_word,
                summary.end_word
            ),
        });
    }

    let mut records = Vec::with_capacity(total);
    for i in 0..total {
        records.push(daf.read_f64_at_word(file_data, data_start_word + i));
    }

    Ok(ChebSegment {
        init,
        intlen,
        rsize,
        ncoeff,
        n_records,
        components,
        records,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic SPK Type 2 segment buffer with a single record
    /// holding constant coefficients per axis.
    fn build_const_type2_segment(x: f64, y: f64, z: f64) -> (Vec<u8>, Daf, Summary) {
        let ncoeff = 1;
        let rsize = 2 + 3 * ncoeff;
        let n_records = 1;
        let n_words = rsize * n_records + 4; // + trailer (init, intlen, rsize, n)
        let mut buf = vec![0u8; n_words * 8];
        let write = |word: usize, v: f64, b: &mut [u8]| {
            b[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
        };
        // Record words 1..=5
        write(1, 0.5, &mut buf); // mid
        write(2, 0.5, &mut buf); // radius
        write(3, x, &mut buf);
        write(4, y, &mut buf);
        write(5, z, &mut buf);
        // Trailer at words 6..=9
        write(6, 0.0, &mut buf); // init
        write(7, 1.0, &mut buf); // intlen
        write(8, rsize as f64, &mut buf);
        write(9, n_records as f64, &mut buf);

        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let summary = Summary {
            start_et: 0.0,
            end_et: 1.0,
            target_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_word: 1,
            end_word: 9,
        };
        (buf, daf, summary)
    }

    #[test]
    fn type2_constant_polynomial_recovers_position() {
        let (buf, daf, summary) = build_const_type2_segment(1.0, 2.0, 3.0);
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        let state = seg.evaluate(0.5).unwrap();
        assert_eq!(state[0], 1.0);
        assert_eq!(state[1], 2.0);
        assert_eq!(state[2], 3.0);
        // Constant polynomial has zero analytic derivative.
        assert_eq!(state[3], 0.0);
        assert_eq!(state[4], 0.0);
        assert_eq!(state[5], 0.0);
    }

    #[test]
    fn type2_linear_polynomial_recovers_velocity() {
        // Build a 1-record segment with ncoeff=2, components=3.
        // p(tau) = c0 + c1 * T1(tau) = c0 + c1 * tau
        // Position at tau is c0 + c1*tau, velocity is c1/radius.
        let ncoeff = 2;
        let rsize = 2 + 3 * ncoeff; // 8
        let n_records = 1;
        let n_words = rsize * n_records + 4; // 12
        let mut buf = vec![0u8; n_words * 8];
        let write = |word: usize, v: f64, b: &mut [u8]| {
            b[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
        };
        let radius = 2.0_f64;
        write(1, 1.0, &mut buf); // mid
        write(2, radius, &mut buf); // radius
                                    // x = 10 + 4*tau ; y = 0 ; z = -1 + 5*tau
        write(3, 10.0, &mut buf);
        write(4, 4.0, &mut buf);
        write(5, 0.0, &mut buf);
        write(6, 0.0, &mut buf);
        write(7, -1.0, &mut buf);
        write(8, 5.0, &mut buf);
        // Trailer at 9..=12
        write(9, 0.0, &mut buf); // init
        write(10, 4.0, &mut buf); // intlen
        write(11, rsize as f64, &mut buf);
        write(12, n_records as f64, &mut buf);

        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let summary = Summary {
            start_et: 0.0,
            end_et: 4.0,
            target_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_word: 1,
            end_word: 12,
        };
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        // Pick et = 1.0 (record midpoint); tau = 0
        let s = seg.evaluate(1.0).unwrap();
        assert!((s[0] - 10.0).abs() < 1e-12);
        assert!((s[1]).abs() < 1e-12);
        assert!((s[2] + 1.0).abs() < 1e-12);
        assert!((s[3] - 4.0 / radius).abs() < 1e-12);
        assert!((s[4]).abs() < 1e-12);
        assert!((s[5] - 5.0 / radius).abs() < 1e-12);
    }

    #[test]
    fn type3_recovers_independent_velocity_polynomial() {
        // Single record, ncoeff=1, components=6. Position constants =
        // (10, 20, 30); velocity constants = (1, 2, 3).
        let ncoeff = 1;
        let rsize = 2 + 6 * ncoeff; // 8
        let n_records = 1;
        let n_words = rsize * n_records + 4;
        let mut buf = vec![0u8; n_words * 8];
        let write = |word: usize, v: f64, b: &mut [u8]| {
            b[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
        };
        write(1, 0.5, &mut buf); // mid
        write(2, 0.5, &mut buf); // radius
        write(3, 10.0, &mut buf);
        write(4, 20.0, &mut buf);
        write(5, 30.0, &mut buf);
        write(6, 1.0, &mut buf);
        write(7, 2.0, &mut buf);
        write(8, 3.0, &mut buf);
        write(9, 0.0, &mut buf); // init
        write(10, 1.0, &mut buf); // intlen
        write(11, rsize as f64, &mut buf);
        write(12, n_records as f64, &mut buf);
        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let summary = Summary {
            start_et: 0.0,
            end_et: 1.0,
            target_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type: 3,
            start_word: 1,
            end_word: 12,
        };
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        let s = seg.evaluate(0.5).unwrap();
        assert_eq!(s, [10.0, 20.0, 30.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn out_of_coverage_returned_for_epoch_outside_range() {
        let (buf, daf, summary) = build_const_type2_segment(1.0, 2.0, 3.0);
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        let err = seg.evaluate(2.0).unwrap_err();
        assert!(matches!(err, SpiceError::OutOfCoverage { .. }));
    }

    #[test]
    fn unsupported_data_type_rejected() {
        let (buf, daf, mut summary) = build_const_type2_segment(1.0, 2.0, 3.0);
        summary.data_type = 13;
        let err = segment_for_summary(&buf, &daf, &summary).unwrap_err();
        assert!(matches!(
            err,
            SpiceError::UnsupportedDataType { data_type: 13 }
        ));
    }

    #[test]
    fn rsize_not_divisible_by_components_rejected() {
        // rsize=6 with components=3 → (6-2)%3=1 → reject
        let n_words = 10;
        let mut buf = vec![0u8; n_words * 8];
        let write = |word: usize, v: f64, b: &mut [u8]| {
            b[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
        };
        write(7, 0.0, &mut buf); // init
        write(8, 1.0, &mut buf); // intlen
        write(9, 6.0, &mut buf); // rsize = 6 (invalid for components=3)
        write(10, 1.0, &mut buf); // n_records
        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let summary = Summary {
            start_et: 0.0,
            end_et: 1.0,
            target_id: 0,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_word: 1,
            end_word: 10,
        };
        let err = segment_for_summary(&buf, &daf, &summary).unwrap_err();
        assert!(matches!(err, SpiceError::Parse { .. }));
    }
}
