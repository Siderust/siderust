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
use crate::time::{JulianDate, TimeContext, TT, UT1};
use crate::qtty::*;

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
///
/// # Panics
///
/// Panics if the bundled `tempoch` ΔT model cannot resolve the supplied
/// epoch. The current model covers approximately 4500 BCE through year
/// 3000 CE; epochs outside that horizon (or non-finite Julian Dates) will
/// trigger an `expect` failure. Domain-restricted callers can use
/// [`tempoch::JulianDate::to_with`] directly to recover a `Result`.
#[inline]
pub fn jd_ut1_from_tt(jd_tt: JulianDate) -> JulianDate {
    let jd_tt: tempoch::JulianDate<crate::time::TT> = jd_tt.into();
    let ut1 = jd_tt
        .to_with::<UT1>(&TimeContext::new())
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    JulianDate::new(ut1.to::<tempoch::JD>().raw().value())
}

/// Compute JD(UT1) from JD(TT) with IERS EOP refinement.
///
/// Routes the conversion through `tempoch`'s context-aware
/// `TT → UT1` chain backed by the bundled IERS `finals2000A.all` series
/// (`TimeContext::with_builtin_eop`). The bundled context applies both the
/// IERS leap-second table (ΔAT) and the daily `dUT1 = UT1 − UTC`
/// interpolation, then the caller's [`EopValues`] is honoured by adding
/// the residual `eop.dut1 − bundled_dut1` correction. In particular:
///
/// * For an [`IersEop`](crate::astro::eop::IersEop) provider that uses the
///   same bundled finals2000A.all data, the residual is ≈ 0 and the
///   result matches tempoch's high-fidelity bundled UT1 path within
///   floating-point precision.
/// * For a [`NullEop`](crate::astro::eop::NullEop) provider (`dut1 = 0`),
///   the residual is `−bundled_dut1`, so the returned value is UTC: the
///   user is asserting `UT1 = UTC` at this epoch, which is the documented
///   semantics of `NullEop`.
/// * For a custom provider with a measured `dut1`, the residual carries
///   only the difference relative to the bundled series — the leap-second
///   chain is still consumed from tempoch's authoritative table.
///
/// This replaces the earlier behaviour of treating `eop` as informational
/// and silently falling back to [`jd_ut1_from_tt`] when `dut1 == 0`: a
/// zero `dut1` is a legitimate EOP assertion, distinct from "no EOP data
/// available, fall back to ΔT". Callers that explicitly want the bundled
/// ΔT model (Stephenson–Houlden, Meeus, IERS-observed 1992–present)
/// should call [`jd_ut1_from_tt`] directly.
///
/// For epochs within the IERS `finals2000A.all` coverage this is accurate
/// to better than 10 ms; outside the observed range tempoch falls back to
/// the same monthly ΔT path used by [`jd_ut1_from_tt`] and the residual
/// from `eop.dut1` is applied on top.
///
/// # Panics
///
/// Panics if `jd_tt` encodes a value that is not finite or falls outside
/// the range supported by the bundled `tempoch` ΔT model (approximately
/// 4500 BCE through year 3000 CE). Callers that need a fallible variant
/// can drive the conversion directly via
/// [`tempoch::JulianDate::to_with`].
#[inline]
pub fn jd_ut1_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> JulianDate {
    let jd_tt_t: tempoch::JulianDate<TT> = jd_tt.into();
    let ctx = TimeContext::with_builtin_eop();
    let ut1 = jd_tt_t
        .to_with::<UT1>(&ctx)
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    let ut1_jd = ut1.to::<tempoch::JD>().raw().value();

    // Honour the caller's `EopValues` by applying the residual
    // `eop.dut1 − bundled_dut1`. We use `ut1_jd` as the UTC-MJD lookup key
    // for the bundled series: |dUT1| ≤ 0.9 s, far below the daily grid
    // spacing of `finals2000A.all`, so the resulting bundle row is the
    // same one tempoch consulted internally.
    let mjd_query = qtty::Day::new(ut1_jd - 2_400_000.5);
    let bundled_dut1 = ctx
        .ut1_minus_utc(mjd_query)
        .map(|s| s.value())
        .unwrap_or(0.0);
    let residual_days = (eop.dut1.value() - bundled_dut1) / 86_400.0;

    JulianDate::new(ut1_jd + residual_days)
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
    fn jd_ut1_eop_zero_dut1_returns_utc() {
        // A zero `dut1` is a legitimate EOP assertion that UT1 = UTC at this
        // epoch. The function must honour the EOP value rather than silently
        // falling back to the ΔT model: the previous behaviour discarded the
        // caller's intent in this case, which is the bug fixed here.
        //
        // `jd_ut1_from_tt_eop(jd, EopValues::default())` therefore returns
        // UTC. We can recover UTC from the bundled UT1 path by subtracting
        // the bundled `ut1_minus_utc`, which is exactly what the function
        // does internally.
        let eop = EopValues::default(); // dut1 = 0 ⇒ UT1 ≡ UTC
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);

        let jd_tt_t: tempoch::JulianDate<TT> = jd().into();
        let ctx = TimeContext::with_builtin_eop();
        let bundled_ut1 = jd_tt_t.to_with::<UT1>(&ctx).expect("bundled UT1");
        let bundled_ut1_jd = bundled_ut1.to::<tempoch::JD>().raw().value();
        let bundled_dut1 = ctx
            .ut1_minus_utc(qtty::Day::new(bundled_ut1_jd - 2_400_000.5))
            .map(|s| s.value())
            .unwrap_or(0.0);
        let expected_utc_jd = bundled_ut1_jd - bundled_dut1 / 86_400.0;

        assert!(
            (jd_eop.value() - expected_utc_jd).abs() < 1e-12,
            "with dut1 = 0, UT1 must equal UTC, got {} vs UTC {}",
            jd_eop.value(),
            expected_utc_jd
        );

        // And the EOP path must NOT collapse to the ΔT-based fallback
        // exposed via `jd_ut1_from_tt`: the two should differ by roughly the
        // bundled dUT1 at this epoch.
        let jd_dt = jd_ut1_from_tt(jd());
        let diff_sec = (jd_eop.value() - jd_dt.value()).abs() * 86_400.0;
        assert!(
            diff_sec > 0.05,
            "EOP(dut1=0) and ΔT models should differ measurably, got {diff_sec}s"
        );
    }

    #[test]
    fn jd_ut1_eop_nonzero_dut1_applies_correction() {
        // For a non-zero user dUT1, the result must equal the bundled UT1
        // adjusted by the residual (user_dut1 − bundled_dut1) so that the
        // overall offset from UTC is exactly the supplied dUT1.
        use crate::qtty::Seconds;
        let dut1_s = 0.3_f64;
        let eop = EopValues {
            dut1: Seconds::new(dut1_s),
            ..Default::default()
        };
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);

        let jd_tt_t: tempoch::JulianDate<TT> = jd().into();
        let ctx = TimeContext::with_builtin_eop();
        let bundled_ut1 = jd_tt_t.to_with::<UT1>(&ctx).expect("bundled UT1");
        let bundled_ut1_jd = bundled_ut1.to::<tempoch::JD>().raw().value();
        let bundled_dut1 = ctx
            .ut1_minus_utc(qtty::Day::new(bundled_ut1_jd - 2_400_000.5))
            .map(|s| s.value())
            .unwrap_or(0.0);
        let expected_utc_jd = bundled_ut1_jd - bundled_dut1 / 86_400.0;
        let recovered_dut1_s = (jd_eop.value() - expected_utc_jd) * 86_400.0;

        // f64 precision at JD ~2.45e6 limits the recoverable resolution to
        // roughly 10 µs; sub-microsecond agreement is unrealistic here.
        assert!(
            (recovered_dut1_s - dut1_s).abs() < 1e-4,
            "UT1-UTC should equal the supplied dUT1 ({dut1_s}s), got {recovered_dut1_s}s"
        );
    }

    #[test]
    fn jd_ut1_eop_negative_dut1() {
        use crate::qtty::Seconds;
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
    fn gmst_from_tt_eop_null_eop_uses_utc_axis() {
        // With NullEop (dut1 = 0), `gmst_from_tt_eop` treats UT1 = UTC, so its
        // value differs from `gmst_from_tt` (ΔT-driven UT1 = true UT1) by the
        // bundled `dUT1` worth of Earth rotation: ω · |bundled_dut1|. Around
        // J2000 the bundled dUT1 is ~−0.3 s ⇒ a few × 10⁻⁵ rad. The exact
        // magnitude depends on the embedded EOP table; we only assert that
        // the two paths *do* differ measurably (the previous behaviour
        // collapsed the EOP path onto the ΔT path).
        let eop = EopValues::default();
        let gmst_eop = gmst_from_tt_eop(jd(), &eop);
        let gmst_dt = gmst_from_tt(jd());

        assert!(
            gmst_eop.value() >= 0.0 && gmst_eop.value() < TAU,
            "EOP-driven GMST out of [0, 2π): {}",
            gmst_eop.value()
        );
        let diff = (gmst_eop.value() - gmst_dt.value()).abs();
        let diff = diff.min((TAU - diff).abs());
        assert!(
            diff > 1e-6,
            "EOP(dut1=0) and ΔT-based GMST should differ by ~ω·dUT1, got {diff} rad"
        );
    }

    #[test]
    fn gmst_from_tt_eop_with_nonzero_dut1() {
        use crate::qtty::Seconds;
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
