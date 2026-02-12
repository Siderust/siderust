// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Segment lookup and coordinate evaluation for DE441 data.
//!
//! Given a Julian Date (TDB), locates the correct Chebyshev record and
//! evaluates the polynomial to produce position `[x, y, z]` in km (ICRF)
//! and optionally velocity `[vx, vy, vz]` in km/day.
//!
//! ## Design — `SegmentDescriptor`
//!
//! All per-body metadata (epoch, interval length, coefficient count, record
//! accessor) is bundled into [`SegmentDescriptor`] so that the evaluation
//! routines are called with a single typed handle instead of 6+ loose
//! parameters.  Static descriptors for each embedded body are exported
//! from the [`data`](super::data) sub-module.

use crate::coordinates::frames::ICRF;
use crate::time::JulianDate;
use affn::{Displacement, Velocity};
use qtty::*;

/// Velocity unit: km/day.
type KmPerDay = Per<Kilometer, Day>;
/// Typed km/day quantity.
type KmPerDayQ = qtty::Quantity<KmPerDay>;

/// Seconds per day — sourced from `qtty` (single source of truth).
const SECONDS_PER_DAY: f64 = qtty::time::SECONDS_PER_DAY;

// ── Segment descriptor ──────────────────────────────────────────────────

/// All metadata needed to evaluate a single DE441 body segment.
///
/// This bundles the per-body constants that were previously passed as
/// six loose parameters (`init`, `intlen`, `ncoeff`, `rsize`,
/// `n_records`, `record_fn`), improving call-site ergonomics and
/// eliminating the risk of mismatched arguments.
pub struct SegmentDescriptor {
    /// Initial epoch of the segment (TDB seconds past J2000).
    pub init: f64,
    /// Length of each Chebyshev sub-interval (seconds).
    pub intlen: f64,
    /// Number of Chebyshev coefficients per coordinate (x, y, z).
    pub ncoeff: usize,
    /// Number of coefficient records in the segment.
    pub n_records: usize,
    /// Accessor for the `i`-th coefficient record.
    pub record_fn: fn(usize) -> &'static [f64],
}

// ── Time conversion ─────────────────────────────────────────────────────

/// Convert a Julian Date (TDB) to TDB seconds past J2000.
#[inline]
fn jd_to_et(jd_tdb: JulianDate) -> Seconds {
    (jd_tdb - JulianDate::J2000).to::<Second>()
}

// ── Internal: locate record + normalise time ────────────────────────────

/// Shared helper: given an epoch return the record, tau, and radius.
#[inline]
fn locate(seg: &SegmentDescriptor, jd_tdb: JulianDate) -> (&'static [f64], f64, f64) {
    let et = jd_to_et(jd_tdb).value();
    let idx = ((et - seg.init) / seg.intlen) as usize;
    let idx = idx.min(seg.n_records - 1);
    let record = (seg.record_fn)(idx);
    let mid = record[0];
    let radius = record[1];
    let tau = (et - mid) / radius;
    (record, tau, radius)
}

/// Extract X/Y/Z coefficient slices from a record.
#[inline]
fn xyz_coeffs(record: &[f64], ncoeff: usize) -> (&[f64], &[f64], &[f64]) {
    (
        &record[2..2 + ncoeff],
        &record[2 + ncoeff..2 + 2 * ncoeff],
        &record[2 + 2 * ncoeff..2 + 3 * ncoeff],
    )
}

// ── Public evaluation API ───────────────────────────────────────────────

impl SegmentDescriptor {
    /// Evaluate position in km (ICRF) at Julian Date (TDB).
    ///
    /// Returns a typed `Displacement<ICRF, Kilometer>` carrying both the
    /// reference frame (ICRF) and the length unit (km) in the type.
    #[inline]
    pub fn position(&self, jd_tdb: JulianDate) -> Displacement<ICRF, Kilometer> {
        let (record, tau, _) = locate(self, jd_tdb);
        let (cx, cy, cz) = xyz_coeffs(record, self.ncoeff);
        Displacement::new(
            Kilometers::new(cheby::evaluate(cx, tau)),
            Kilometers::new(cheby::evaluate(cy, tau)),
            Kilometers::new(cheby::evaluate(cz, tau)),
        )
    }

    /// Evaluate velocity in km/day (ICRF) at Julian Date (TDB).
    ///
    /// Returns a typed `Velocity<ICRF, Per<Kilometer, Day>>`.
    ///
    /// The Chebyshev derivative gives `df/dτ`; we multiply by `dτ/dt = 1/radius`
    /// (km/s) then by `SECONDS_PER_DAY` to obtain km/day.
    #[inline]
    pub fn velocity(&self, jd_tdb: JulianDate) -> Velocity<ICRF, KmPerDay> {
        let (record, tau, radius) = locate(self, jd_tdb);
        let (cx, cy, cz) = xyz_coeffs(record, self.ncoeff);
        let scale = SECONDS_PER_DAY / radius;
        Velocity::new(
            KmPerDayQ::new(cheby::evaluate_derivative(cx, tau) * scale),
            KmPerDayQ::new(cheby::evaluate_derivative(cy, tau) * scale),
            KmPerDayQ::new(cheby::evaluate_derivative(cz, tau) * scale),
        )
    }

    /// Evaluate both position and velocity in one pass.
    ///
    /// Returns `(Displacement<ICRF, Kilometer>, Velocity<ICRF, KmPerDay>)`.
    #[inline]
    pub fn position_velocity(
        &self,
        jd_tdb: JulianDate,
    ) -> (Displacement<ICRF, Kilometer>, Velocity<ICRF, KmPerDay>) {
        let (record, tau, radius) = locate(self, jd_tdb);
        let (cx, cy, cz) = xyz_coeffs(record, self.ncoeff);
        let scale = SECONDS_PER_DAY / radius;

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
}
