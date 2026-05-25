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
//! This module covers the SPK data types relevant to the JPL DE planetary
//! kernels and common spacecraft kernels exercised by POD:
//!
//! * **Type 2** — Chebyshev polynomials for position only; velocity is
//!   recovered analytically from the position polynomial derivative.
//! * **Type 3** — Chebyshev polynomials for position **and** velocity,
//!   with independent coefficient sets.
//! * **Type 9** — unequally-spaced states interpolated with Lagrange
//!   polynomials.
//! * **Type 13** — unequally-spaced states interpolated with Hermite
//!   polynomials using stored velocities as first derivatives.
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

/// A typed SPK segment.
#[derive(Debug, Clone)]
pub enum SpkSegment {
    /// SPK Type 2 — Chebyshev polynomial for position; velocity is
    /// derived analytically.
    Type2(ChebSegment),
    /// SPK Type 3 — Chebyshev polynomials for position and velocity
    /// stored side by side (no analytic differentiation needed).
    Type3(ChebSegment),
    /// SPK Type 9 — unequally-spaced state records with Lagrange
    /// interpolation.
    Type9(DiscreteStatesSegment),
    /// SPK Type 13 — unequally-spaced state records with Hermite
    /// interpolation.
    Type13(DiscreteStatesSegment),
}

impl SpkSegment {
    /// Coverage start (TDB seconds past J2000) of the segment.
    pub fn start_tdb_seconds(&self) -> f64 {
        match self {
            SpkSegment::Type2(c) | SpkSegment::Type3(c) => c.init,
            SpkSegment::Type9(s) | SpkSegment::Type13(s) => s.start_tdb_seconds(),
        }
    }

    /// Coverage end (TDB seconds past J2000) of the segment.
    pub fn end_tdb_seconds(&self) -> f64 {
        match self {
            SpkSegment::Type2(c) | SpkSegment::Type3(c) => c.init + c.intlen * c.n_records as f64,
            SpkSegment::Type9(s) | SpkSegment::Type13(s) => s.end_tdb_seconds(),
        }
    }

    /// Borrow the underlying Chebyshev segment data.
    pub fn cheb(&self) -> &ChebSegment {
        match self {
            SpkSegment::Type2(c) | SpkSegment::Type3(c) => c,
            SpkSegment::Type9(_) | SpkSegment::Type13(_) => {
                panic!("SPK segment is not a Chebyshev Type 2/3 segment")
            }
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
            SpkSegment::Type9(s) => s.evaluate_lagrange(et),
            SpkSegment::Type13(s) => s.evaluate_hermite(et),
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

/// Unequally-spaced state records shared by SPK Type 9 and Type 13.
#[derive(Debug, Clone)]
pub struct DiscreteStatesSegment {
    /// SPK data type (`9` or `13`).
    pub data_type: i32,
    /// Lagrange degree for Type 9, or Hermite sample count for Type 13.
    pub interpolation_size: usize,
    /// Number of state records.
    pub n_records: usize,
    /// Flat `[x, y, z, vx, vy, vz]` records in km and km/s.
    pub states: Vec<f64>,
    /// Epoch for each state record, in TDB seconds past J2000.
    pub epochs: Vec<f64>,
    /// Optional SPK epoch directory entries.
    pub epoch_directory: Vec<f64>,
}

impl DiscreteStatesSegment {
    /// Coverage start (TDB seconds past J2000).
    pub fn start_tdb_seconds(&self) -> f64 {
        self.epochs.first().copied().unwrap_or(f64::NAN)
    }

    /// Coverage end (TDB seconds past J2000).
    pub fn end_tdb_seconds(&self) -> f64 {
        self.epochs.last().copied().unwrap_or(f64::NAN)
    }

    fn state_record(&self, idx: usize) -> Result<[f64; 6], SpiceError> {
        let start = idx.checked_mul(6).ok_or_else(|| SpiceError::Parse {
            message: format!("state index overflow for record {idx}"),
        })?;
        let slice = self
            .states
            .get(start..start + 6)
            .ok_or_else(|| SpiceError::Parse {
                message: format!("state record {idx} outside Type {} data", self.data_type),
            })?;
        Ok([slice[0], slice[1], slice[2], slice[3], slice[4], slice[5]])
    }

    fn locate_or_window(&self, et: f64, count: usize) -> Result<WindowSelection, SpiceError> {
        if !et.is_finite() {
            return Err(SpiceError::Corrupted {
                message: format!("epoch is not finite: {et}"),
            });
        }
        let start = self.start_tdb_seconds();
        let end = self.end_tdb_seconds();
        if self.n_records == 0 || et < start || et > end {
            return Err(SpiceError::OutOfCoverage {
                target: 0,
                center: 0,
                epoch_tdb_seconds: et,
                start_tdb_seconds: start,
                end_tdb_seconds: end,
            });
        }
        match self.epochs.binary_search_by(|epoch| {
            epoch
                .partial_cmp(&et)
                .unwrap_or(core::cmp::Ordering::Less)
        }) {
            Ok(idx) => Ok(WindowSelection::Exact(idx)),
            Err(insertion) => {
                if count == 0 || count > self.n_records {
                    return Err(SpiceError::Parse {
                        message: format!(
                            "Type {} interpolation window {count} invalid for {} records",
                            self.data_type, self.n_records
                        ),
                    });
                }
                let left = count / 2;
                let mut first = insertion.saturating_sub(left);
                if first + count > self.n_records {
                    first = self.n_records - count;
                }
                Ok(WindowSelection::Window { first, count })
            }
        }
    }

    fn evaluate_lagrange(&self, et: f64) -> Result<[f64; 6], SpiceError> {
        let count = self.interpolation_size + 1;
        match self.locate_or_window(et, count)? {
            WindowSelection::Exact(idx) => self.state_record(idx),
            WindowSelection::Window { first, count } => {
                let epochs = &self.epochs[first..first + count];
                let mut out = [0.0_f64; 6];
                for (component, output) in out.iter_mut().enumerate() {
                    let mut values = Vec::with_capacity(count);
                    for idx in first..first + count {
                        values.push(self.states[idx * 6 + component]);
                    }
                    *output = lagrange_value(epochs, &values, et)?;
                }
                Ok(out)
            }
        }
    }

    fn evaluate_hermite(&self, et: f64) -> Result<[f64; 6], SpiceError> {
        let count = self.interpolation_size;
        match self.locate_or_window(et, count)? {
            WindowSelection::Exact(idx) => self.state_record(idx),
            WindowSelection::Window { first, count } => {
                let epochs = &self.epochs[first..first + count];
                let mut out = [0.0_f64; 6];
                for axis in 0..3 {
                    let mut positions = Vec::with_capacity(count);
                    let mut velocities = Vec::with_capacity(count);
                    for idx in first..first + count {
                        positions.push(self.states[idx * 6 + axis]);
                        velocities.push(self.states[idx * 6 + axis + 3]);
                    }
                    let (position, velocity) = hermite_value_and_derivative(
                        epochs,
                        &positions,
                        &velocities,
                        et,
                    )?;
                    out[axis] = position;
                    out[axis + 3] = velocity;
                }
                Ok(out)
            }
        }
    }
}

enum WindowSelection {
    Exact(usize),
    Window { first: usize, count: usize },
}

/// Read an SPK segment summary into a typed [`SpkSegment`].
///
/// Supports SPK Type 2, Type 3, Type 9, and Type 13; any other type is rejected with
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
        9 => Ok(SpkSegment::Type9(read_discrete_states(
            file_data, daf, summary, 9,
        )?)),
        13 => Ok(SpkSegment::Type13(read_discrete_states(
            file_data, daf, summary, 13,
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

fn read_discrete_states(
    file_data: &[u8],
    daf: &Daf,
    summary: &Summary,
    data_type: i32,
) -> Result<DiscreteStatesSegment, SpiceError> {
    let n_words = summary
        .end_word
        .checked_sub(summary.start_word)
        .and_then(|delta| delta.checked_add(1))
        .ok_or_else(|| SpiceError::Parse {
            message: format!(
                "invalid SPK Type {data_type} word range [{}..{}]",
                summary.start_word, summary.end_word
            ),
        })?;
    if n_words < 9 {
        return Err(SpiceError::Parse {
            message: format!("SPK Type {data_type} segment too short: {n_words} words"),
        });
    }
    let mut words = Vec::with_capacity(n_words);
    for i in 0..n_words {
        words.push(daf.read_f64_at_word(file_data, summary.start_word + i));
    }

    let n_records_f = words[n_words - 1];
    let descriptor_f = words[n_words - 2];
    if !n_records_f.is_finite() || n_records_f < 1.0 {
        return Err(SpiceError::Parse {
            message: format!("SPK Type {data_type} invalid record count {n_records_f}"),
        });
    }
    if !descriptor_f.is_finite() || descriptor_f < 0.0 {
        return Err(SpiceError::Parse {
            message: format!("SPK Type {data_type} invalid interpolation descriptor {descriptor_f}"),
        });
    }
    let n_records = n_records_f as usize;
    let descriptor = descriptor_f as usize;
    let state_len = n_records.checked_mul(6).ok_or_else(|| SpiceError::Parse {
        message: format!("SPK Type {data_type} state count overflow: {n_records}"),
    })?;
    let epoch_start = state_len;
    let epoch_end = epoch_start + n_records;
    if epoch_end + 2 > words.len() {
        return Err(SpiceError::Parse {
            message: format!(
                "SPK Type {data_type} payload too short for {n_records} records"
            ),
        });
    }
    let states = words[..state_len].to_vec();
    let epochs = words[epoch_start..epoch_end].to_vec();
    let epoch_directory = words[epoch_end..words.len() - 2].to_vec();
    if states
        .iter()
        .chain(epochs.iter())
        .chain(epoch_directory.iter())
        .any(|v| !v.is_finite())
    {
        return Err(SpiceError::Corrupted {
            message: format!("SPK Type {data_type} contains non-finite state or epoch data"),
        });
    }
    if !epochs.windows(2).all(|w| w[0] < w[1]) {
        return Err(SpiceError::Parse {
            message: format!("SPK Type {data_type} epochs are not strictly increasing"),
        });
    }

    let interpolation_size = match data_type {
        9 => {
            let degree = descriptor;
            if degree == 0 || degree + 1 > n_records || degree > 31 {
                return Err(SpiceError::Parse {
                    message: format!(
                        "SPK Type 9 invalid Lagrange degree {degree} for {n_records} records"
                    ),
                });
            }
            degree
        }
        13 => {
            let samples = descriptor + 1;
            if samples < 2 || samples > n_records || samples > 16 {
                return Err(SpiceError::Parse {
                    message: format!(
                        "SPK Type 13 invalid Hermite sample count {samples} for {n_records} records"
                    ),
                });
            }
            samples
        }
        _ => unreachable!(),
    };

    Ok(DiscreteStatesSegment {
        data_type,
        interpolation_size,
        n_records,
        states,
        epochs,
        epoch_directory,
    })
}

fn lagrange_value(xs: &[f64], ys: &[f64], x: f64) -> Result<f64, SpiceError> {
    if xs.len() != ys.len() || xs.is_empty() {
        return Err(SpiceError::Parse {
            message: "Lagrange interpolation requires matching non-empty inputs".into(),
        });
    }
    let mut value = 0.0_f64;
    for i in 0..xs.len() {
        let mut basis = 1.0_f64;
        for j in 0..xs.len() {
            if i == j {
                continue;
            }
            let denom = xs[i] - xs[j];
            if denom.abs() <= f64::EPSILON {
                return Err(SpiceError::Corrupted {
                    message: "duplicate epochs in Lagrange interpolation window".into(),
                });
            }
            basis *= (x - xs[j]) / denom;
        }
        value += ys[i] * basis;
    }
    Ok(value)
}

fn hermite_value_and_derivative(
    xs: &[f64],
    ys: &[f64],
    dydx: &[f64],
    x: f64,
) -> Result<(f64, f64), SpiceError> {
    if xs.len() != ys.len() || xs.len() != dydx.len() || xs.is_empty() {
        return Err(SpiceError::Parse {
            message: "Hermite interpolation requires matching non-empty inputs".into(),
        });
    }
    let n = xs.len();
    let order = 2 * n;
    let mut z = vec![0.0_f64; order];
    let mut table = vec![vec![0.0_f64; order]; order];

    for i in 0..n {
        let even = 2 * i;
        let odd = even + 1;
        z[even] = xs[i];
        z[odd] = xs[i];
        table[even][0] = ys[i];
        table[odd][0] = ys[i];
        table[odd][1] = dydx[i];
        if i == 0 {
            table[even][1] = dydx[i];
        } else {
            let denom = z[even] - z[even - 1];
            if denom.abs() <= f64::EPSILON {
                return Err(SpiceError::Corrupted {
                    message: "duplicate epochs in Hermite interpolation window".into(),
                });
            }
            table[even][1] = (table[even][0] - table[even - 1][0]) / denom;
        }
    }

    for i in 2..order {
        for j in 2..=i {
            let denom = z[i] - z[i - j];
            if denom.abs() <= f64::EPSILON {
                return Err(SpiceError::Corrupted {
                    message: "duplicate epochs in Hermite divided differences".into(),
                });
            }
            table[i][j] = (table[i][j - 1] - table[i - 1][j - 1]) / denom;
        }
    }

    let coeffs: Vec<f64> = (0..order).map(|i| table[i][i]).collect();
    let mut value = coeffs[order - 1];
    let mut derivative = 0.0_f64;
    for i in (0..order - 1).rev() {
        derivative = derivative * (x - z[i]) + value;
        value = value * (x - z[i]) + coeffs[i];
    }
    Ok((value, derivative))
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
        summary.data_type = 21;
        let err = segment_for_summary(&buf, &daf, &summary).unwrap_err();
        assert!(matches!(
            err,
            SpiceError::UnsupportedDataType { data_type: 21 }
        ));
    }

    fn build_discrete_segment(data_type: i32, descriptor: f64) -> (Vec<u8>, Daf, Summary) {
        let n_records = 4usize;
        let n_words = n_records * 6 + n_records + 2;
        let mut buf = vec![0u8; n_words * 8];
        let write = |word: usize, v: f64, b: &mut [u8]| {
            b[(word - 1) * 8..word * 8].copy_from_slice(&v.to_le_bytes());
        };
        for i in 0..n_records {
            let t = i as f64;
            let base = i * 6 + 1;
            // Linear motion: x=t, y=2t, z=3t; velocity is constant.
            write(base, t, &mut buf);
            write(base + 1, 2.0 * t, &mut buf);
            write(base + 2, 3.0 * t, &mut buf);
            write(base + 3, 1.0, &mut buf);
            write(base + 4, 2.0, &mut buf);
            write(base + 5, 3.0, &mut buf);
        }
        let epoch_start = n_records * 6 + 1;
        for i in 0..n_records {
            write(epoch_start + i, i as f64, &mut buf);
        }
        write(n_words - 1, descriptor, &mut buf);
        write(n_words, n_records as f64, &mut buf);
        let daf = Daf {
            nd: 2,
            ni: 6,
            summaries: vec![],
        };
        let summary = Summary {
            start_et: 0.0,
            end_et: 3.0,
            target_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type,
            start_word: 1,
            end_word: n_words,
        };
        (buf, daf, summary)
    }

    #[test]
    fn type9_lagrange_interpolates_state_records() {
        let (buf, daf, summary) = build_discrete_segment(9, 1.0);
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        let s = seg.evaluate(1.5).unwrap();
        assert!((s[0] - 1.5).abs() < 1e-12);
        assert!((s[1] - 3.0).abs() < 1e-12);
        assert!((s[2] - 4.5).abs() < 1e-12);
        assert!((s[3] - 1.0).abs() < 1e-12);
        assert!((s[4] - 2.0).abs() < 1e-12);
        assert!((s[5] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn type13_hermite_interpolates_position_and_velocity() {
        // Type 13 stores sample count minus one in the descriptor.
        let (buf, daf, summary) = build_discrete_segment(13, 1.0);
        let seg = segment_for_summary(&buf, &daf, &summary).unwrap();
        let s = seg.evaluate(1.5).unwrap();
        assert!((s[0] - 1.5).abs() < 1e-12);
        assert!((s[1] - 3.0).abs() < 1e-12);
        assert!((s[2] - 4.5).abs() < 1e-12);
        assert!((s[3] - 1.0).abs() < 1e-12);
        assert!((s[4] - 2.0).abs() < 1e-12);
        assert!((s[5] - 3.0).abs() < 1e-12);
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
