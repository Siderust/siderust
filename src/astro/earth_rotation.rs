// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Rotation Helpers
//!
//! Utilities for computing UT1 Julian Dates from TT, and properly
//! separated GMST with distinct UT1 and TT arguments.
//!
//! These helpers centralise the TT→UT1 conversion logic so that all
//! call sites (topocentric, horizontal, stellar) use a consistent,
//! IAU-compliant approach.

use crate::astro::eop::EopValues;
use crate::astro::sidereal::gmst_iau2006;
use crate::time::{JulianDate, UT};
use qtty::*;

/// Compute JD(UT1) from JD(TT) using tempoch's ΔT model.
///
/// This uses the piecewise ΔT model in `tempoch` (Stephenson & Houlden,
/// Meeus table, IERS-observed 1992–2025, linear extrapolation),
/// accurate to ~0.5 s for modern epochs.
///
/// The result is stored as a `JulianDate` (which has TT semantics in
/// tempoch's type system), but its numeric value is on the UT1 axis.
/// This matches the SOFA convention where `iauGmst06` takes UT1 as a
/// plain JD double.
#[inline]
pub fn jd_ut1_from_tt(jd_tt: JulianDate) -> JulianDate {
    // tempoch's UT scale: from_jd_tt subtracts ΔT via 3-iteration fixed point.
    // .value() gives us the raw JD number on the UT1 axis.
    let ut1_value = jd_tt.to::<UT>().value();
    JulianDate::new(ut1_value)
}

/// Compute JD(UT1) from JD(TT) with IERS EOP refinement.
///
/// Uses the IERS leap-second table (via `tempoch::tai_minus_utc`) and the
/// EOP `dUT1 = UT1 − UTC` value to compute UT1 via the exact chain:
///
/// ```text
/// UTC = TT − (ΔAT + 32.184 s)
/// UT1 = UTC + dUT1
/// ```
///
/// When `dUT1` is zero (e.g. from [`NullEop`](crate::astro::eop::NullEop)),
/// this falls back to the tempoch ΔT model via [`jd_ut1_from_tt`].
///
/// For epochs within the IERS `finals2000A.all` table, this is accurate
/// to ~1 ms (limited by the EOP table's sampling and interpolation).
#[inline]
pub fn jd_ut1_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> JulianDate {
    // If EOP dUT1 is exactly zero, fall back to the ΔT model.
    // This handles NullEop gracefully and avoids a subtle bias:
    // the leap-second chain with dUT1=0 gives UTC, not UT1.
    if eop.dut1.value() == 0.0 {
        return jd_ut1_from_tt(jd_tt);
    }

    // Exact chain: TT → UTC → UT1
    //   TT = TAI + 32.184 s
    //   TAI = UTC + ΔAT            (ΔAT = cumulative leap seconds)
    //   ⇒ UTC = TT − (ΔAT + 32.184) / 86400        (in days)
    //   UT1 = UTC + dUT1 / 86400                     (in days)
    //
    // We approximate JD(UTC) to look up ΔAT, then refine.
    const TT_MINUS_TAI: f64 = 32.184; // seconds
    let approx_utc = jd_tt.value() - (37.0 + TT_MINUS_TAI) / 86_400.0;
    let dat = tempoch::tai_minus_utc(approx_utc);
    let jd_utc = jd_tt.value() - (dat + TT_MINUS_TAI) / 86_400.0;
    let jd_ut1 = jd_utc + eop.dut1.value() / 86_400.0;

    JulianDate::new(jd_ut1)
}

/// Compute GMST with proper UT1/TT separation.
///
/// This is a convenience wrapper around [`gmst_iau2006`] that automatically
/// computes JD(UT1) from JD(TT) using tempoch's ΔT model.
///
/// # Arguments
///
/// * `jd_tt` - Julian Date on the TT time scale
///
/// # Returns
///
/// Greenwich Mean Sidereal Time in radians
#[inline]
pub fn gmst_from_tt(jd_tt: JulianDate) -> Radians {
    let jd_ut1 = jd_ut1_from_tt(jd_tt);
    gmst_iau2006(jd_ut1, jd_tt)
}

/// Compute GMST with proper UT1/TT separation and EOP refinement.
///
/// Uses the IERS EOP data to refine the UT1 computation.
#[inline]
pub fn gmst_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> Radians {
    let jd_ut1 = jd_ut1_from_tt_eop(jd_tt, eop);
    gmst_iau2006(jd_ut1, jd_tt)
}

/// Compute GMST using the default embedded IERS EOP table.
///
/// This is the recommended function for general use. It provides
/// IAU 2006-compliant GMST with proper UT1/TT separation using
/// tempoch's ΔT model.
#[inline]
pub fn gmst_default(jd_tt: JulianDate) -> Radians {
    gmst_from_tt(jd_tt)
}
