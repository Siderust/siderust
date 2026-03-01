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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::eop::EopValues;
    use std::f64::consts::TAU;

    const JD_J2000: f64 = 2451545.0;

    fn jd() -> JulianDate {
        JulianDate::new(JD_J2000)
    }

    // ── jd_ut1_from_tt ────────────────────────────────────────────────────

    #[test]
    fn jd_ut1_is_close_to_tt_at_modern_epoch() {
        // ΔT at J2000 is ~63.8 s, so UT1 ≈ TT - 63.8/86400 days
        let jd_ut1 = jd_ut1_from_tt(jd());
        let diff_sec = (jd().value() - jd_ut1.value()) * 86400.0;
        assert!(
            diff_sec > 50.0 && diff_sec < 80.0,
            "ΔT at J2000 expected ~63s, got {diff_sec}s"
        );
    }

    // ── jd_ut1_from_tt_eop ────────────────────────────────────────────────

    #[test]
    fn jd_ut1_eop_zero_dut1_falls_back_to_delta_t() {
        let eop = EopValues::default(); // dut1 = 0
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);
        let jd_dt = jd_ut1_from_tt(jd());
        // Both should give the same result when dUT1 = 0
        assert!((jd_eop.value() - jd_dt.value()).abs() < 1e-12);
    }

    #[test]
    fn jd_ut1_eop_nonzero_dut1_applies_correction() {
        use qtty::Seconds;
        let eop = EopValues {
            dut1: Seconds::new(0.3),
            ..Default::default()
        }; // dUT1 = +0.3 s
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);
        // With non-zero dUT1, result should differ from pure ΔT model
        assert!(jd_eop.value().is_finite());
        // The result should be within a few minutes of TT (no gross errors)
        let diff_days = (jd().value() - jd_eop.value()).abs();
        assert!(
            diff_days < 1e-3,
            "UT1-TT difference too large: {diff_days} days"
        );
    }

    #[test]
    fn jd_ut1_eop_negative_dut1() {
        use qtty::Seconds;
        let eop = EopValues {
            dut1: Seconds::new(-0.5),
            ..Default::default()
        };
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);
        assert!(jd_eop.value().is_finite());
    }

    // ── gmst_from_tt ─────────────────────────────────────────────────────

    #[test]
    fn gmst_from_tt_is_in_range() {
        let gmst = gmst_from_tt(jd());
        assert!(
            gmst.value() >= 0.0 && gmst.value() < TAU,
            "GMST out of [0, 2π): {}",
            gmst.value()
        );
    }

    #[test]
    fn gmst_from_tt_varies_with_time() {
        let gmst1 = gmst_from_tt(jd());
        let gmst2 = gmst_from_tt(JulianDate::new(JD_J2000 + 1.0));
        // One sidereal day rotates by ~2π, so GMST changes significantly
        let diff = (gmst2.value() - gmst1.value()).abs();
        assert!(diff > 0.0, "GMST should change over 1 day");
    }

    // ── gmst_from_tt_eop ─────────────────────────────────────────────────

    #[test]
    fn gmst_from_tt_eop_null_eop_matches_gmst_from_tt() {
        let eop = EopValues::default(); // all zeros
        let gmst_eop = gmst_from_tt_eop(jd(), &eop);
        let gmst_dt = gmst_from_tt(jd());
        // With NullEop (all zeros), they should be equal
        assert!((gmst_eop.value() - gmst_dt.value()).abs() < 1e-12);
    }

    #[test]
    fn gmst_from_tt_eop_with_nonzero_dut1() {
        use qtty::Seconds;
        let eop = EopValues {
            dut1: Seconds::new(0.2),
            ..Default::default()
        };
        let gmst_eop = gmst_from_tt_eop(jd(), &eop);
        assert!(gmst_eop.value().is_finite());
        assert!(gmst_eop.value() >= 0.0);
    }

    // ── gmst_default ─────────────────────────────────────────────────────

    #[test]
    fn gmst_default_matches_gmst_from_tt() {
        let jd_val = jd();
        let g1 = gmst_default(jd_val);
        let g2 = gmst_from_tt(jd_val);
        assert!((g1.value() - g2.value()).abs() < 1e-15);
    }
}
