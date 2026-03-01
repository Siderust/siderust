// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared segment lookup and coordinate evaluation for DE4xx datasets.
//!
//! This module provides two descriptor types:
//! - [`SegmentDescriptor`]: compile-time data via function pointer to `&'static [f64]`.
//! - [`DynSegmentDescriptor`]: runtime data backed by a heap-allocated `Vec<f64>`.
//!
//! Both share the same Chebyshev evaluation logic.

use crate::coordinates::frames::ICRF;
use crate::time::JulianDate;
use affn::{Displacement, Velocity};
use qtty::*;

type KmPerDay = Per<Kilometer, Day>;
type KmPerDayQ = qtty::Quantity<KmPerDay>;

const SECONDS_PER_DAY: f64 = qtty::time::SECONDS_PER_DAY;

// ═══════════════════════════════════════════════════════════════════════════
// Shared helpers (used by both SegmentDescriptor and DynSegmentDescriptor)
// ═══════════════════════════════════════════════════════════════════════════

#[inline]
fn jd_to_et(jd_tdb: JulianDate) -> Seconds {
    (jd_tdb - JulianDate::J2000).to::<Second>()
}

#[inline]
fn xyz_coeffs(record: &[f64], ncoeff: usize) -> (&[f64], &[f64], &[f64]) {
    (
        &record[2..2 + ncoeff],
        &record[2 + ncoeff..2 + 2 * ncoeff],
        &record[2 + 2 * ncoeff..2 + 3 * ncoeff],
    )
}

/// Compute the record index and Chebyshev parameter tau from segment metadata.
#[inline]
fn locate_params(
    init: Seconds,
    intlen: Seconds,
    n_records: usize,
    jd_tdb: JulianDate,
) -> (usize, f64, Seconds) {
    let et = jd_to_et(jd_tdb);
    let idx = ((et - init) / intlen).value() as usize;
    let idx = idx.min(n_records - 1);
    (idx, et.value(), intlen)
}

/// Given a record slice, compute (tau, radius).
#[inline]
fn record_tau(record: &[f64], et: f64) -> (f64, Seconds) {
    let mid = Seconds::new(record[0]);
    let radius = Seconds::new(record[1]);
    let tau = (et - mid.value()) / radius.value();
    (tau, radius)
}

/// Evaluate position from a record slice.
#[inline]
fn eval_position(record: &[f64], ncoeff: usize, tau: f64) -> Displacement<ICRF, Kilometer> {
    let (cx, cy, cz) = xyz_coeffs(record, ncoeff);
    Displacement::new(
        Kilometers::new(cheby::evaluate(cx, tau)),
        Kilometers::new(cheby::evaluate(cy, tau)),
        Kilometers::new(cheby::evaluate(cz, tau)),
    )
}

/// Evaluate velocity from a record slice.
#[inline]
fn eval_velocity(
    record: &[f64],
    ncoeff: usize,
    tau: f64,
    radius: Seconds,
) -> Velocity<ICRF, KmPerDay> {
    let (cx, cy, cz) = xyz_coeffs(record, ncoeff);
    let scale = SECONDS_PER_DAY / radius.value();
    Velocity::new(
        KmPerDayQ::new(cheby::evaluate_derivative(cx, tau) * scale),
        KmPerDayQ::new(cheby::evaluate_derivative(cy, tau) * scale),
        KmPerDayQ::new(cheby::evaluate_derivative(cz, tau) * scale),
    )
}

/// Evaluate both position and velocity from a record slice.
#[inline]
fn eval_both(
    record: &[f64],
    ncoeff: usize,
    tau: f64,
    radius: Seconds,
) -> (Displacement<ICRF, Kilometer>, Velocity<ICRF, KmPerDay>) {
    let (cx, cy, cz) = xyz_coeffs(record, ncoeff);
    let scale = SECONDS_PER_DAY / radius.value();

    let (px, vx) = cheby::evaluate_both(cx, tau);
    let (py, vy) = cheby::evaluate_both(cy, tau);
    let (pz, vz) = cheby::evaluate_both(cz, tau);

    (
        Displacement::new(
            Kilometers::new(px),
            Kilometers::new(py),
            Kilometers::new(pz),
        ),
        Velocity::new(
            KmPerDayQ::new(vx * scale),
            KmPerDayQ::new(vy * scale),
            KmPerDayQ::new(vz * scale),
        ),
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// SegmentDescriptor — compile-time embedded data (via function pointer)
// ═══════════════════════════════════════════════════════════════════════════

/// All metadata needed to evaluate one body segment (compile-time data).
///
/// The `record_fn` field points to a function that returns a `&'static [f64]`
/// slice into data embedded at compile time via `include_bytes!`.
pub struct SegmentDescriptor {
    /// Initial epoch of the segment (TDB seconds past J2000).
    pub init: Seconds,
    /// Length of each Chebyshev sub-interval (seconds).
    pub intlen: Seconds,
    /// Number of Chebyshev coefficients per coordinate (x, y, z).
    pub ncoeff: usize,
    /// Number of coefficient records in the segment.
    pub n_records: usize,
    /// Accessor for the `i`-th coefficient record.
    pub record_fn: fn(usize) -> &'static [f64],
}

#[inline]
fn locate(seg: &SegmentDescriptor, jd_tdb: JulianDate) -> (&'static [f64], f64, Seconds) {
    let (idx, et, intlen) = locate_params(seg.init, seg.intlen, seg.n_records, jd_tdb);
    let record = (seg.record_fn)(idx);
    let (tau, radius) = record_tau(record, et);
    let _ = intlen;
    (record, tau, radius)
}

impl SegmentDescriptor {
    /// Evaluate position in km (ICRF) at Julian Date (TDB).
    #[inline]
    pub fn position(&self, jd_tdb: JulianDate) -> Displacement<ICRF, Kilometer> {
        let (record, tau, _) = locate(self, jd_tdb);
        eval_position(record, self.ncoeff, tau)
    }

    /// Evaluate velocity in km/day (ICRF) at Julian Date (TDB).
    #[inline]
    pub fn velocity(&self, jd_tdb: JulianDate) -> Velocity<ICRF, KmPerDay> {
        let (record, tau, radius) = locate(self, jd_tdb);
        eval_velocity(record, self.ncoeff, tau, radius)
    }

    /// Evaluate both position and velocity in one pass.
    #[inline]
    pub fn position_velocity(
        &self,
        jd_tdb: JulianDate,
    ) -> (Displacement<ICRF, Kilometer>, Velocity<ICRF, KmPerDay>) {
        let (record, tau, radius) = locate(self, jd_tdb);
        eval_both(record, self.ncoeff, tau, radius)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// DynSegmentDescriptor — runtime heap-backed data
// ═══════════════════════════════════════════════════════════════════════════

/// Segment descriptor for runtime-loaded data (heap-backed).
///
/// Unlike [`SegmentDescriptor`], the coefficient data is stored in a `Vec<f64>`
/// owned by this struct (or shared via `Arc`). This enables loading BSP files
/// at runtime without embedding them at compile time.
pub struct DynSegmentDescriptor {
    /// Initial epoch of the segment (TDB seconds past J2000).
    pub init: Seconds,
    /// Length of each Chebyshev sub-interval (seconds).
    pub intlen: Seconds,
    /// Number of Chebyshev coefficients per coordinate (x, y, z).
    pub ncoeff: usize,
    /// Doubles per record (2 + 3 × ncoeff).
    pub rsize: usize,
    /// Number of coefficient records in the segment.
    pub n_records: usize,
    /// Flattened coefficient data: n_records × rsize f64 values.
    pub data: Vec<f64>,
}

impl DynSegmentDescriptor {
    /// Construct from parsed SPK segment data.
    pub fn from_spk(seg: &crate::data::spk::SegmentData) -> Self {
        Self {
            init: Seconds::new(seg.init),
            intlen: Seconds::new(seg.intlen),
            ncoeff: seg.ncoeff,
            rsize: seg.rsize,
            n_records: seg.n_records,
            data: seg.records.clone(),
        }
    }

    /// Get the `i`-th record as a slice.
    #[inline]
    fn record(&self, i: usize) -> &[f64] {
        let start = i * self.rsize;
        &self.data[start..start + self.rsize]
    }

    /// Locate the record for a given JD and compute tau.
    #[inline]
    fn locate(&self, jd_tdb: JulianDate) -> (&[f64], f64, Seconds) {
        let (idx, et, _) = locate_params(self.init, self.intlen, self.n_records, jd_tdb);
        let record = self.record(idx);
        let (tau, radius) = record_tau(record, et);
        (record, tau, radius)
    }

    /// Evaluate position in km (ICRF) at Julian Date (TDB).
    #[inline]
    pub fn position(&self, jd_tdb: JulianDate) -> Displacement<ICRF, Kilometer> {
        let (record, tau, _) = self.locate(jd_tdb);
        eval_position(record, self.ncoeff, tau)
    }

    /// Evaluate velocity in km/day (ICRF) at Julian Date (TDB).
    #[inline]
    pub fn velocity(&self, jd_tdb: JulianDate) -> Velocity<ICRF, KmPerDay> {
        let (record, tau, radius) = self.locate(jd_tdb);
        eval_velocity(record, self.ncoeff, tau, radius)
    }

    /// Evaluate both position and velocity in one pass.
    #[inline]
    pub fn position_velocity(
        &self,
        jd_tdb: JulianDate,
    ) -> (Displacement<ICRF, Kilometer>, Velocity<ICRF, KmPerDay>) {
        let (record, tau, radius) = self.locate(jd_tdb);
        eval_both(record, self.ncoeff, tau, radius)
    }
}
