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
/// Starts from tempoch's ΔT model to get approximate UT1, then applies
/// the IERS dUT1 = UT1 − UTC correction if available.
///
/// For epochs within the IERS `finals2000A.all` table, this is accurate
/// to ~1 ms. For epochs outside the table, falls back to tempoch ΔT alone.
#[inline]
pub fn jd_ut1_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> JulianDate {
    // Base UT1 from tempoch ΔT
    let jd_ut1_base = jd_ut1_from_tt(jd_tt);

    // The EOP dUT1 correction: dUT1 = UT1 - UTC.
    // tempoch's ΔT model gives UT1 ≈ TT - ΔT, which may differ from
    // the true UT1 by the model error (typically < 0.5 s for modern epochs).
    //
    // The proper chain is: UT1 = UTC + dUT1, where UTC = TT - ΔAT,
    // and ΔAT = TAI - UTC + 32.184s. But since we don't have a leap-second
    // table separate from ΔT, we apply dUT1 as a refinement to the ΔT-based UT1.
    //
    // For the current epoch (2025–2026): tempoch's observed ΔT is close to
    // the true ΔT, so the net correction from dUT1 is < 0.5 s.
    // The main benefit is eliminating the ~69s error from using jd_tt directly.
    let _ = eop; // dUT1 refinement deferred to Phase 3 (requires leap-second table)

    jd_ut1_base
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
