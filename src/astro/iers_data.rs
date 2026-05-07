// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IERS Earth Orientation Parameter Lookup
//!
//! ## Scientific scope
//!
//! The **International Earth Rotation and Reference Systems Service** (IERS)
//! publishes daily measurements and short-range predictions of Earth Orientation
//! Parameters (EOP): polar motion `xp`, `yp`; the difference between Universal
//! Time and Coordinated Universal Time `ΔUT1 = UT1 − UTC`; and celestial-pole
//! offsets `dX`, `dY` in the IAU 2000A/2006 precession-nutation model.  These
//! quantities are required for transforming between terrestrial and celestial
//! reference frames at sub-arcsecond accuracy.
//!
//! This compatibility module wraps the historical auto-generated
//! `iers_eop_data.rs` table and provides binary-search interpolation keyed on
//! Modified Julian Date in UTC. New Siderust code should use
//! [`crate::astro::eop::IersEop`], which consumes the active/bundled EOP data
//! owned by `tempoch`.
//!
//! ## Technical scope
//!
//! - [`EopEntry`] — one row of the IERS table: `mjd`, `xp`, `yp`, `dut1`,
//!   `dx`, `dy` (all `f64` raw values from the generated code).
//! - [`EOP_TABLE`] — static slice of `EopEntry` rows.
//! - [`lookup`] — binary-search + linear interpolation keyed on [`Days`] (MJD
//!   in UTC).  Returns `None` when the epoch is outside the table's coverage.
//!
//! The `mjd` argument is typed as [`Days`] (`Quantity<Day>`) so callers cannot
//! accidentally pass a Julian Date (JD) or Julian Centuries value.  The
//! internal comparison against `EopEntry.mjd` (which remains `f64`) extracts
//! `.value()` at the function boundary.
//!
//! ## References
//!
//! - Petit, G., & Luzum, B. (Eds.) (2010). *IERS Conventions (2010)*. IERS
//!   Technical Note 36, §5.1. Verlag des BKG, Frankfurt.
//! - IERS EOP Data Center. <https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html>

#[allow(clippy::approx_constant)]
#[rustfmt::skip]
mod iers_eop_data {
    include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/generated/iers_eop_data.rs"));
}

pub use iers_eop_data::{EopEntry, EOP_TABLE};

use crate::qtty::Days;

/// Look up EOP values by MJD (UTC), interpolating linearly between
/// daily entries when the MJD falls between two table rows.
///
/// # Arguments
///
/// - `table`: slice of [`EopEntry`] rows, typically [`EOP_TABLE`].
/// - `mjd`: Modified Julian Date in UTC scale, typed as [`Days`].  Callers
///   must subtract the JD epoch offset themselves:
///   `Days::new(jd_utc - 2_400_000.5)`.
///
/// # Returns
///
/// * `Some(entry)` when `mjd` is within the closed interval
///   `[first.mjd, last.mjd]` covered by `table`.
/// * `None` in two cases:
///     - The supplied `table` slice is empty (i.e. the embedded EOP table
///       has not been generated, or a custom table was passed in empty).
///     - `mjd` falls strictly outside the `[first, last]` interval. The
///       caller should return or propagate a missing-EOP error unless a
///       documented `NullEop` approximation was selected explicitly.
///       `iers_data::lookup` deliberately does **not** extrapolate, since
///       IERS entries beyond the prediction horizon are unreliable and would
///       silently degrade UT1/polar-motion accuracy without a clear signal.
///
/// Within the covered interval the function performs daily linear
/// interpolation and never returns `None` — `None` is therefore a
/// reliable indicator that the requested epoch is out of range.
///
/// # Examples
///
/// ```
/// use siderust::astro::iers_data::{EOP_TABLE, lookup};
/// use siderust::qtty::Days;
///
/// // MJD 51544.5 ≈ J2000.0 in UTC
/// let result = lookup(&EOP_TABLE, Days::new(51544.5));
/// // Table may or may not contain this epoch depending on the build-time data.
/// let _ = result;
/// ```
pub fn lookup(table: &[EopEntry], mjd: Days) -> Option<EopEntry> {
    let mjd = mjd.value();
    if table.is_empty() {
        return None;
    }

    let first = table.first().unwrap().mjd;
    let last = table.last().unwrap().mjd;

    if mjd < first || mjd > last {
        return None;
    }

    // Binary search for the interval containing `mjd`.
    let idx = table.partition_point(|e| e.mjd < mjd);

    if idx == 0 {
        return Some(table[0]);
    }
    if idx >= table.len() {
        return Some(*table.last().unwrap());
    }

    let lo = &table[idx - 1];
    let hi = &table[idx];

    // Exact match?
    if (hi.mjd - mjd).abs() < 1e-10 {
        return Some(*hi);
    }
    if (lo.mjd - mjd).abs() < 1e-10 {
        return Some(*lo);
    }

    // Linear interpolation
    let dt = hi.mjd - lo.mjd;
    let t = (mjd - lo.mjd) / dt;

    Some(EopEntry {
        mjd,
        xp: lo.xp + t * (hi.xp - lo.xp),
        yp: lo.yp + t * (hi.yp - lo.yp),
        dut1: lo.dut1 + t * (hi.dut1 - lo.dut1),
        dx: lo.dx + t * (hi.dx - lo.dx),
        dy: lo.dy + t * (hi.dy - lo.dy),
    })
}
