// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Rotation Helpers
//!
//! Utilities for computing UT1/UTC Julian Dates from TT and for evaluating GMST
//! with properly separated UT1 and TT arguments.
//!
//! ## Scientific scope
//!
//! Many Earth-rotation quantities (ERA, GMST, GAST, sidereal time) are
//! defined on the **UT1** time scale, while the precession, nutation and
//! polynomial parts are defined on the dynamical **TT** scale. Mixing the
//! two introduces tens-of-seconds-level errors. The conversion `TT → UT1`
//! requires the time-dependent quantity `ΔT = TT − UT1`, which combines the
//! tabulated leap-second history (ΔAT = TAI − UTC) with the IERS-published
//! `dUT1 = UT1 − UTC` series.
//!
//! ## Technical scope
//!
//! **Time-scale ownership:** `tempoch` owns leap-second history, ΔT tables,
//! EOP bundles, and the active time-data status. Siderust provides GMST/sidereal
//! adapters and coordinate-pipeline convenience functions on top.
//!
//! | Function | EOP required? | Notes |
//! |---|---|---|
//! | [`jd_ut1_from_tt`] | no | monthly ΔT model via [`TimeContext::new`] |
//! | [`try_jd_utc_from_tt`] | **yes** | IERS-indexed UTC for EOP lookup keys |
//! | [`jd_utc_from_tt_delta_t`] | no | ΔT + leap-second approximation |
//! | [`jd_ut1_from_tt_eop`] | preferred | refines UT1 when runtime EOP is active |
//! | [`gmst_from_tt`] | no | ΔT-derived UT1 |
//! | [`gmst_default`] | no | alias of [`gmst_from_tt`] (ΔT approximation) |
//! | [`try_gmst_with_eop`] | **yes** | IERS EOP–refined GMST |
//!
//! EOP data is available only when a `tempoch` runtime bundle is active
//! (loaded via `tempoch`'s runtime-data mechanisms). Without runtime EOP,
//! [`try_jd_utc_from_tt`] and [`try_gmst_with_eop`] return
//! [`EopError::NoData`]; they never fabricate UTC by equating UT1 with UTC.
//!
//! The `JulianDate` returned by the UT1/UTC helpers carries the numeric JD
//! value of the requested timescale but is typed as TT for interoperability
//! with the rest of the sidereal-time call sites (matching the SOFA convention
//! for `iauGmst06`).
//!
//! ## References
//!
//! * IERS Conventions (2010), §5.4 (ERA and GAST), §5.5 (GMST)
//! * SOFA routine `iauGmst06`
//! * Stephenson & Houlden (1986); Meeus, *Astronomical Algorithms*, Ch. 10

use crate::astro::eop::{EopError, EopProvider, EopValues, IersEop};
use crate::astro::sidereal::gmst_iau2006;
use crate::qtty::*;
use crate::time::{JulianDate, TimeContext, TT, UT1};

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
    let jd_tt: tempoch::JulianDate<crate::time::TT> = jd_tt;
    let ut1 = jd_tt
        .to_with::<UT1>(&TimeContext::new())
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    crate::time::JulianDate::new(ut1.to::<tempoch::JD>().raw().value())
}

/// Compute JD(UTC) from JD(TT) using the active IERS EOP bundle.
///
/// This is the UTC axis required for IERS EOP table lookups. The conversion
/// uses a fixed-point iteration on `UT1 − UTC` from the runtime-loaded
/// `tempoch` EOP series.
///
/// # Errors
///
/// Returns [`EopError::NoData`] when no runtime EOP bundle covers the epoch.
/// Callers that only need a ΔT-based civil-time approximation should use
/// [`jd_utc_from_tt_delta_t`].
#[inline]
pub fn try_jd_utc_from_tt(jd_tt: JulianDate) -> Result<JulianDate, EopError> {
    let jd_tt_t: tempoch::JulianDate<TT> = jd_tt;
    let ctx = TimeContext::with_builtin_eop();
    let ut1 = jd_tt_t
        .to_with::<UT1>(&ctx)
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    let ut1_jd = ut1.to::<tempoch::JD>().raw().value();
    let mut utc_jd = ut1_jd;

    for _ in 0..3 {
        let mjd_utc = Days::new(utc_jd - 2_400_000.5);
        let Some(eop) = tempoch::eop::builtin_eop_at(mjd_utc) else {
            return Err(EopError::NoData {
                jd_utc: utc_jd,
                mjd_utc: mjd_utc.value(),
            });
        };
        utc_jd = ut1_jd - eop.ut1_minus_utc.to::<Day>().value();
    }

    Ok(crate::time::JulianDate::new(utc_jd))
}

/// Compute JD(UTC) from JD(TT) without requiring runtime EOP.
///
/// Uses the ΔT model to obtain JD(UT1) and treats that axis as UTC, neglecting
/// the sub-second `UT1 − UTC` offset when no daily EOP series is available
/// (typically < 0.9 s). This is adequate for coarse civil-time indexing but
/// not for IERS EOP table keys — use [`try_jd_utc_from_tt`] there.
#[inline]
pub fn jd_utc_from_tt_delta_t(jd_tt: JulianDate) -> JulianDate {
    jd_ut1_from_tt(jd_tt)
}

/// Compute JD(UTC) from JD(TT).
///
/// # Deprecated behaviour
///
/// This is a compatibility wrapper around [`jd_utc_from_tt_delta_t`]. It does
/// **not** consult the runtime IERS EOP bundle and therefore must not be used
/// as an EOP table index when sub-second UT1 accuracy is required. For IERS
/// lookups use [`try_jd_utc_from_tt`].
#[inline]
pub fn jd_utc_from_tt(jd_tt: JulianDate) -> JulianDate {
    jd_utc_from_tt_delta_t(jd_tt)
}

/// Compute JD(UT1) from JD(TT) with IERS EOP refinement.
///
/// Routes the conversion through `tempoch`'s context-aware `TT → UT1` chain.
/// When a runtime EOP bundle is active, the bundled daily `dUT1 = UT1 − UTC`
/// interpolation is applied first; the caller's [`EopValues`] is then honoured
/// by adding the residual `eop.dut1 − bundled_dut1`.
///
/// Without runtime EOP the bundled `dUT1` term is unavailable and the
/// residual reduces to the caller-supplied `eop.dut1` relative to zero.
///
/// # Panics
///
/// Panics if `jd_tt` encodes a value that is not finite or falls outside
/// the range supported by the bundled `tempoch` ΔT model (approximately
/// 4500 BCE through year 3000 CE).
#[inline]
pub fn jd_ut1_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> JulianDate {
    let jd_tt_t: tempoch::JulianDate<TT> = jd_tt;
    let ctx = TimeContext::with_builtin_eop();
    let ut1 = jd_tt_t
        .to_with::<UT1>(&ctx)
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    let ut1_jd = ut1.to::<tempoch::JD>().raw().value();

    let mjd_query = qtty::Day::new(ut1_jd - 2_400_000.5);
    let bundled_dut1 = ctx
        .ut1_minus_utc(mjd_query)
        .map(|s| s.value())
        .unwrap_or(0.0);
    let residual_days = (eop.dut1.value() - bundled_dut1) / 86_400.0;

    crate::time::JulianDate::new(ut1_jd + residual_days)
}

/// Compute GMST with proper UT1/TT separation.
///
/// Convenience wrapper around [`gmst_iau2006`] that computes JD(UT1) from
/// JD(TT) using tempoch's ΔT model.
#[inline]
pub fn gmst_from_tt(jd_tt: JulianDate) -> Radians {
    let jd_ut1 = jd_ut1_from_tt(jd_tt);
    gmst_iau2006(jd_ut1, jd_tt)
}

/// Compute GMST with proper UT1/TT separation and EOP refinement.
#[inline]
pub fn gmst_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> Radians {
    let jd_ut1 = jd_ut1_from_tt_eop(jd_tt, eop);
    gmst_iau2006(jd_ut1, jd_tt)
}

/// Compute GMST using the ΔT model (no runtime EOP required).
///
/// This is the recommended default when a runtime IERS EOP bundle is not
/// loaded. Accuracy is limited by the monthly ΔT table (~0.5 s in UT1 for
/// modern epochs). For IERS-refined GMST use [`try_gmst_with_eop`].
#[inline]
pub fn gmst_default(jd_tt: JulianDate) -> Radians {
    gmst_from_tt(jd_tt)
}

/// Compute GMST with IERS EOP refinement when a runtime bundle is active.
///
/// # Errors
///
/// Returns [`EopError::NoData`] when no runtime EOP bundle covers the epoch.
#[inline]
pub fn try_gmst_with_eop(jd_tt: JulianDate) -> Result<Radians, EopError> {
    let jd_utc = try_jd_utc_from_tt(jd_tt)?;
    let eop = IersEop::new().try_eop_at(jd_utc)?;
    Ok(gmst_from_tt_eop(jd_tt, &eop))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::eop::EopValues;
    use std::f64::consts::TAU;

    const JD_J2000: f64 = tempoch::J2000_JD_TT_DAY.value();

    fn jd() -> JulianDate {
        crate::time::JulianDate::new(JD_J2000)
    }

    // ── jd_ut1_from_tt ────────────────────────────────────────────────────

    #[test]
    fn jd_ut1_is_close_to_tt_at_modern_epoch() {
        let jd_ut1 = jd_ut1_from_tt(jd());
        let diff_sec = (jd().raw().value() - jd_ut1.raw().value()) * 86400.0;
        assert!(
            diff_sec > 50.0 && diff_sec < 80.0,
            "ΔT at J2000 expected ~63s, got {diff_sec}s"
        );
    }

    // ── try_jd_utc_from_tt ────────────────────────────────────────────────

    #[test]
    fn try_jd_utc_returns_no_data_without_runtime_eop() {
        assert!(matches!(
            try_jd_utc_from_tt(jd()),
            Err(EopError::NoData { .. })
        ));
    }

    #[test]
    fn jd_utc_delta_t_differs_from_tt_at_j2000() {
        let jd_utc = jd_utc_from_tt_delta_t(jd());
        let diff_sec = (jd().raw().value() - jd_utc.raw().value()) * 86_400.0;
        assert!(
            diff_sec > 60.0 && diff_sec < 70.0,
            "TT-UTC at J2000 should include leap seconds + 32.184s, got {diff_sec}s"
        );
    }

    // ── jd_ut1_from_tt_eop ────────────────────────────────────────────────

    #[test]
    fn jd_ut1_eop_zero_dut1_returns_utc() {
        use crate::qtty::Seconds;
        let eop_zero = EopValues::default();
        let eop_nonzero = EopValues {
            dut1: Seconds::new(0.5),
            ..Default::default()
        };
        let jd_zero = jd_ut1_from_tt_eop(jd(), &eop_zero);
        let jd_nonzero = jd_ut1_from_tt_eop(jd(), &eop_nonzero);
        let diff_sec = (jd_nonzero.raw().value() - jd_zero.raw().value()).abs() * 86_400.0;
        assert!(
            (diff_sec - 0.5).abs() < 0.01,
            "nonzero dut1 should shift UT1 by ~0.5s, got {diff_sec}s"
        );
    }

    #[test]
    fn jd_ut1_eop_nonzero_dut1_applies_correction() {
        use crate::qtty::Seconds;
        let dut1_s = 0.3_f64;
        let eop = EopValues {
            dut1: Seconds::new(dut1_s),
            ..Default::default()
        };
        let jd_eop = jd_ut1_from_tt_eop(jd(), &eop);
        let jd_delta_t = jd_ut1_from_tt(jd());
        let recovered_dut1_s = (jd_eop.raw().value() - jd_delta_t.raw().value()) * 86_400.0;
        assert!(
            (recovered_dut1_s - dut1_s).abs() < 1e-4,
            "residual dut1 should equal supplied value ({dut1_s}s), got {recovered_dut1_s}s"
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
        assert!(jd_eop.raw().value().is_finite());
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
        let gmst2 = gmst_from_tt(crate::time::JulianDate::new(JD_J2000 + 1.0));
        let diff = (gmst2.value() - gmst1.value()).abs();
        assert!(diff > 0.0, "GMST should change over 1 day");
    }

    // ── gmst_from_tt_eop ─────────────────────────────────────────────────

    #[test]
    fn gmst_from_tt_eop_null_eop_uses_utc_axis() {
        use crate::qtty::Seconds;
        let eop_zero = EopValues::default();
        let eop_nonzero = EopValues {
            dut1: Seconds::new(0.3),
            ..Default::default()
        };
        let gmst_zero = gmst_from_tt_eop(jd(), &eop_zero);
        let gmst_nonzero = gmst_from_tt_eop(jd(), &eop_nonzero);

        assert!(
            gmst_zero.value() >= 0.0 && gmst_zero.value() < TAU,
            "EOP-driven GMST out of [0, 2π): {}",
            gmst_zero.value()
        );
        let diff = (gmst_nonzero.value() - gmst_zero.value()).abs();
        let diff = diff.min((TAU - diff).abs());
        assert!(
            diff > 1e-6,
            "EOP dut1=0.3s should shift GMST by ~2×10⁻⁵ rad, got {diff} rad"
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
        assert!(gmst_eop.is_finite());
        assert!(gmst_eop >= Radians::new(0.0));
    }

    // ── gmst_default / try_gmst_with_eop ──────────────────────────────────

    #[test]
    fn gmst_default_matches_delta_t_path() {
        let jd_val = jd();
        assert!((gmst_default(jd_val).value() - gmst_from_tt(jd_val).value()).abs() < 1e-15);
    }

    #[test]
    fn try_gmst_with_eop_requires_runtime_bundle() {
        assert!(matches!(
            try_gmst_with_eop(jd()),
            Err(EopError::NoData { .. })
        ));
    }
}
