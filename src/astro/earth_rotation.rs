// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Rotation Helpers
//!
//! Utilities for computing UT1 Julian Dates from TT and for evaluating GMST
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
//! | [`jd_ut1_from_tt_eop`] | preferred | refines UT1 when runtime EOP is active |
//! | [`gmst_from_tt`] | no | ΔT-derived UT1 |
//! | [`try_gmst_with_eop`] | **yes** | IERS EOP–refined GMST |
//!
//! EOP data is available only when a `tempoch` runtime bundle is active
//! (loaded via `tempoch`'s runtime-data mechanisms). Without runtime EOP,
//! [`try_gmst_with_eop`] returns
//! [`EopError::NoData`]. UTC coordinates come from tempoch's native conversion.
//!
//! UT1 arguments and results carry `JulianDate<UT1>`; dynamical arguments
//! carry `JulianDate<TT>`.
//!
//! ## References
//!
//! * IERS Conventions (2010), §5.4 (ERA and GAST), §5.5 (GMST)
//! * SOFA routine `iauGmst06`
//! * Stephenson & Houlden (1986); Meeus, *Astronomical Algorithms*, Ch. 10

use crate::astro::eop::{EopError, EopProvider, EopValues, IersEop};
use crate::astro::sidereal::gmst_iau2006;
use crate::qtty::*;
use crate::time::{JulianDate, TimeContext, UT1};

/// Compute JD(UT1) from JD(TT) using tempoch's ΔT model.
///
/// This uses the piecewise ΔT model in `tempoch` (Stephenson & Houlden,
/// Meeus table, IERS-observed 1992–2025, linear extrapolation),
/// accurate to ~0.5 s for modern epochs.
///
/// # Panics
///
/// Panics if the bundled `tempoch` ΔT model cannot resolve the supplied
/// epoch. The current model covers approximately 4500 BCE through year
/// 3000 CE; epochs outside that horizon (or non-finite Julian Dates) will
/// trigger an `expect` failure. Domain-restricted callers can use
/// [`tempoch::JulianDate::to_with`] directly to recover a `Result`.
#[inline]
pub fn jd_ut1_from_tt(jd_tt: JulianDate) -> tempoch::JulianDate<UT1> {
    let ut1 = jd_tt
        .to_with::<UT1>(&TimeContext::new())
        .expect("TT->UT1 conversion should succeed within the bundled model horizon");
    ut1.to::<tempoch::JD>()
}

/// Compute JD(UT1) from JD(TT) with IERS EOP refinement.
///
/// The supplied DUT1 defines UT1 relative to native civil UTC. The context
/// supplies TAI-UTC, including at leap-second instants. No monthly ΔT
/// estimate is added to the supplied EOP.
///
/// # Panics
///
/// Panics if the UTC offset or resulting UT1 coordinate cannot be represented.
#[inline]
pub fn jd_ut1_from_tt_eop(jd_tt: JulianDate, eop: &EopValues) -> tempoch::JulianDate<UT1> {
    // A supplied DUT1 defines UT1 relative to civil UTC, independently of
    // the monthly ΔT model. Use its atomic coordinate and the context-owned
    // offset without passing a leap instant through a scalar UTC JD.
    let ctx = TimeContext::new().allow_pre_definition_utc();
    let utc = jd_tt.to::<tempoch::UTC>();
    let offset = ctx
        .tai_minus_utc(utc)
        .expect("TT epoch must have a representable civil UTC date");
    let (hi, lo) = utc.to::<tempoch::TAI>().raw_seconds_pair();
    tempoch::Time::<UT1>::try_from_raw_j2000_seconds_split(hi, lo - offset + eop.dut1)
        .expect("UT1 coordinate must be finite")
        .to::<tempoch::JD>()
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

/// Compute GMST with IERS EOP refinement when a runtime bundle is active.
///
/// # Errors
///
/// Returns [`EopError::NoData`] when no runtime EOP bundle covers the epoch.
#[inline]
pub fn try_gmst_with_eop(jd_tt: JulianDate) -> Result<Radians, EopError> {
    let jd_utc = jd_tt
        .to::<tempoch::UTC>()
        .to_with::<tempoch::JD>(&TimeContext::new())?;
    let eop = IersEop::new().try_eop_at(jd_utc)?;
    Ok(gmst_from_tt_eop(jd_tt, &eop))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::eop::EopValues;
    use core::f64::consts::TAU;

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
        let jd_utc = jd()
            .to::<tempoch::UTC>()
            .to_with::<tempoch::JD>(&TimeContext::new())
            .unwrap();
        let recovered_dut1_s = (jd_eop.raw().value()
            - jd_utc.try_raw_with(&TimeContext::new()).unwrap().value())
            * 86_400.0;
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

    // ── try_gmst_with_eop ──────────────────────────────────

    #[test]
    fn try_gmst_with_eop_requires_runtime_bundle() {
        assert!(matches!(
            try_gmst_with_eop(jd()),
            Err(EopError::NoData { .. })
        ));
    }
}
