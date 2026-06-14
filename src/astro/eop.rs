// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Orientation Parameters (EOP)
//!
//! Loading, representation and provider abstraction for Earth Orientation
//! Parameters used throughout the Earth-rotation pipeline.
//!
//! ## Scientific scope
//!
//! Earth Orientation Parameters connect the celestial (GCRS) and terrestrial
//! (ITRS) reference frames by accounting for irregularities in Earth's
//! rotation that cannot be predicted from theory alone. The relevant
//! quantities are:
//!
//! | Parameter | Description |
//! |-----------|-------------|
//! | **UT1 − UTC** | Difference between Universal Time and Coordinated Universal Time |
//! | **xₚ, yₚ** | Pole coordinates, position of the CIP in the ITRS |
//! | **dX, dY** | Celestial pole offsets, corrections to the IAU precession-nutation model |
//! | **LOD** | Length Of Day excess (deviation from 86400 SI seconds) |
//!
//! These are published by the IERS in Bulletins A/B and the
//! `finals2000A.all` series, and are required at the sub-arcsecond level for
//! VLBI, GNSS and astrometric reductions.
//!
//! ## Technical scope
//!
//! The [`EopProvider`] trait abstracts over data sources (zero-correction
//! placeholder, bundled IERS tables, custom feeds). [`EopValues`] holds the
//! parameters in their natural IERS units (`Seconds`, `Arcseconds`,
//! `MilliArcseconds`); converting to radians for computation goes through
//! the standard `qtty` machinery.
//!
//! **Time-scale contract:** all Julian Dates passed to EOP lookups are in
//! **UTC**. IERS Bulletin A/B tables are indexed by UTC civil date; passing
//! TT or UT1 would introduce ~70 s of error in the interpolation index and
//! corrupt the returned `dUT1`, `xp`, `yp` values.
//!
//! ## References
//!
//! * IERS Conventions (2010), Chapter 5
//! * IERS Technical Note 36

use crate::qtty::*;
use crate::time::JulianDate;

/// A set of Earth Orientation Parameters at a given epoch.
///
/// Field units match IERS publications for natural construction from
/// Bulletin A/B data:
/// - `dut1` in [`Seconds`] (UT1 − UTC)
/// - `xp`, `yp` in [`Arcseconds`] (pole coordinates)
/// - `dx`, `dy` in [`crate::qtty::MilliArcseconds`] (celestial pole offsets)
///
/// Use `.to::<Radian>()` on any angular field to get radians for computation.
#[derive(Debug, Clone, Copy)]
pub struct EopValues {
    /// UT1 − UTC.
    pub dut1: Seconds,
    /// Pole x-coordinate.
    pub xp: Arcseconds,
    /// Pole y-coordinate.
    pub yp: Arcseconds,
    /// Celestial pole offset dX.
    pub dx: MilliArcseconds,
    /// Celestial pole offset dY.
    pub dy: MilliArcseconds,
    /// Length-of-day excess (deviation from 86400 SI seconds).
    pub lod: Seconds,
}

impl Default for EopValues {
    fn default() -> Self {
        Self {
            dut1: Seconds::new(0.0),
            xp: Arcseconds::new(0.0),
            yp: Arcseconds::new(0.0),
            dx: MilliArcseconds::new(0.0),
            dy: MilliArcseconds::new(0.0),
            lod: Seconds::new(0.0),
        }
    }
}

impl EopValues {
    /// Convert UTC Julian Date to UT1 Julian Date using this EOP's dUT1.
    #[inline]
    pub fn jd_ut1(&self, jd_utc: JulianDate) -> JulianDate {
        crate::time::JulianDate::new(jd_utc.raw().value() + self.dut1.to::<Day>().value())
    }
}

/// EOP lookup error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EopError {
    /// No EOP entry covers the requested UTC epoch.
    NoData {
        /// Requested UTC Julian Date.
        jd_utc: f64,
        /// Requested UTC Modified Julian Date.
        mjd_utc: f64,
    },
}

impl core::fmt::Display for EopError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::NoData { jd_utc, mjd_utc } => write!(
                f,
                "no Earth Orientation Parameters cover UTC JD {jd_utc} (MJD {mjd_utc})"
            ),
        }
    }
}

impl std::error::Error for EopError {}

/// Trait for providing Earth Orientation Parameters at a given epoch.
///
/// Implementors might read from IERS data files, interpolate tables,
/// or return zero corrections for simplified calculations.
///
/// # Time-scale contract
///
/// **All `jd_*` arguments throughout this trait and its implementations are
/// UTC Julian Dates.**  Callers must not pass TT, TDB, or UT1 Julian Dates
/// without first converting to UTC (via the appropriate Δt tables).
///
/// The reason is that IERS Bulletin A/B tables are indexed by UTC civil date;
/// using the wrong time scale introduces a bias of up to ~70 s (≈ accumulated
/// leap seconds + TT–TAI) in the interpolation, which can shift the looked-up
/// `dUT1` by several milliseconds and shift `xp`/`yp` by a perceptible amount.
pub trait EopProvider {
    /// Fallibly look up EOP values for the given **UTC** Julian Date.
    ///
    /// The default implementation preserves compatibility for providers that
    /// are intentionally total, such as [`NullEop`].
    #[inline]
    fn try_eop_at(&self, jd_utc: JulianDate) -> Result<EopValues, EopError> {
        Ok(self.eop_at(jd_utc))
    }

    /// Look up (or interpolate) EOP values for the given **UTC** Julian Date.
    ///
    /// # Time-scale contract
    /// `jd_utc` **must** be a UTC Julian Date.  Passing a TT or UT1 Julian
    /// Date will silently produce incorrect EOP values due to the offset
    /// between UTC and those scales (up to ~70 s accumulated).
    fn eop_at(&self, jd_utc: JulianDate) -> EopValues;
}

/// Null EOP provider: returns zero for all parameters.
///
/// This is suitable when:
/// - EOP data is unavailable
/// - Accuracy requirements are < 1″ (polar motion effect)
/// - UTC ≈ UT1 is acceptable (|dUT1| < 0.9 s)
#[derive(Debug, Clone, Copy, Default)]
pub struct NullEop;

impl EopProvider for NullEop {
    #[inline]
    fn try_eop_at(&self, _jd_utc: JulianDate) -> Result<EopValues, EopError> {
        Ok(EopValues::default())
    }

    #[inline]
    fn eop_at(&self, _jd_utc: JulianDate) -> EopValues {
        EopValues::default()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// IERS EOP provider (runtime-loaded via tempoch)
// ═══════════════════════════════════════════════════════════════════════════

/// EOP provider backed by tempoch's active IERS `finals2000A.all` bundle.
///
/// Uses the time-data bundle that `tempoch` has loaded at runtime.  EOP data
/// is **not** compiled into the binary; it is loaded on demand via
/// `tempoch`'s runtime-data mechanisms.  Runtime refreshes performed through
/// `tempoch` are immediately visible to newly constructed providers because
/// `IersEop` always reads from tempoch's active bundle.
///
/// The provider interpolates linearly between daily entries for any epoch
/// within the bundle's MJD range. For epochs outside the range,
/// [`IersEop::try_eop_at`] returns [`EopError::NoData`]; the provider never
/// fabricates zero EOP values outside coverage.
#[derive(Debug, Clone, Copy, Default)]
pub struct IersEop;

impl IersEop {
    /// Create an `IersEop` provider.
    ///
    /// EOP data is not embedded at compile time; the provider reads from
    /// whichever `tempoch` time-data bundle is active at the time of lookup.
    #[inline]
    pub fn new() -> Self {
        Self
    }

    /// MJD range covered by the active table.
    ///
    /// # Returns
    ///
    /// * `Some((first_mjd, last_mjd))` when the table has at least one
    ///   entry; the typed [`Days`] values are the first and last MJD entries,
    ///   inclusive on both ends.
    /// * `None` when no runtime EOP bundle has been loaded.
    pub fn mjd_range(&self) -> Option<(Days, Days)> {
        let start = tempoch::eop::eop_start()?;
        let end = tempoch::eop::eop_end()?;
        Some((start, end))
    }

    /// Number of entries in the active table.
    #[inline]
    pub fn len(&self) -> usize {
        match self.mjd_range() {
            Some((start, end)) => (end - start).value().floor() as usize + 1,
            None => 0,
        }
    }

    /// Whether the table is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.mjd_range().is_none()
    }
}

impl EopProvider for IersEop {
    fn try_eop_at(&self, jd_utc: JulianDate) -> Result<EopValues, EopError> {
        let mjd = Days::new(jd_utc.raw().value() - 2_400_000.5);

        let e = tempoch::eop::builtin_eop_at(mjd).ok_or(EopError::NoData {
            jd_utc: jd_utc.raw().value(),
            mjd_utc: mjd.value(),
        })?;

        Ok(EopValues {
            dut1: e.ut1_minus_utc,
            xp: Arcseconds::new(e.pm_xp.map(|v| v.value()).unwrap_or(0.0)),
            yp: Arcseconds::new(e.pm_yp.map(|v| v.value()).unwrap_or(0.0)),
            dx: MilliArcseconds::new(e.dx.map(|v| v.value()).unwrap_or(0.0)),
            dy: MilliArcseconds::new(e.dy.map(|v| v.value()).unwrap_or(0.0)),
            lod: Seconds::new(e.lod.map(|v| v.value() / 1_000.0).unwrap_or(0.0)),
        })
    }

    fn eop_at(&self, jd_utc: JulianDate) -> EopValues {
        self.try_eop_at(jd_utc).unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_eop_returns_zeros() {
        let eop = NullEop;
        let vals = eop.eop_at(crate::J2000);
        assert_eq!(vals.dut1.value(), 0.0);
        assert_eq!(vals.xp.value(), 0.0);
        assert_eq!(vals.yp.value(), 0.0);
        assert_eq!(vals.dx.value(), 0.0);
        assert_eq!(vals.dy.value(), 0.0);
    }

    #[test]
    fn eop_jd_ut1_conversion() {
        let vals = EopValues {
            dut1: Seconds::new(0.35),
            ..Default::default()
        };
        let jd_utc = crate::J2000;
        let jd_ut1 = vals.jd_ut1(jd_utc);
        let diff_days = (jd_ut1.raw() - jd_utc.raw()).value();
        let diff_s = diff_days * 86400.0;
        assert!(
            (diff_s - 0.35).abs() < 1e-3,
            "UT1-UTC = {}s, expected 0.35s",
            diff_s
        );
    }

    #[test]
    fn eop_unit_conversions() {
        let as2rad = std::f64::consts::PI / (180.0 * 3600.0);
        let vals = EopValues {
            xp: Arcseconds::new(0.1),
            yp: Arcseconds::new(0.2),
            dx: MilliArcseconds::new(0.3),
            dy: MilliArcseconds::new(0.4),
            ..Default::default()
        };
        assert!((vals.xp.to::<Radian>().value() - 0.1 * as2rad).abs() < 1e-20);
        assert!((vals.yp.to::<Radian>().value() - 0.2 * as2rad).abs() < 1e-20);
        assert!((vals.dx.to::<Radian>().value() - 0.3e-3 * as2rad).abs() < 1e-20);
        assert!((vals.dy.to::<Radian>().value() - 0.4e-3 * as2rad).abs() < 1e-20);
    }

    #[test]
    fn iers_eop_returns_no_data_without_runtime_bundle() {
        // tempoch no longer embeds EOP at compile time; IersEop.try_eop_at
        // returns Err(NoData) when no runtime bundle has been loaded.
        let eop = IersEop::new();
        let result = eop.try_eop_at(crate::J2000);
        assert!(
            matches!(result, Err(EopError::NoData { .. })),
            "expected NoData without runtime EOP bundle, got {result:?}"
        );
    }

    #[test]
    fn iers_eop_falls_back_to_zeros_without_runtime_bundle() {
        // eop_at() uses the graceful fallback (zero corrections) when no
        // runtime EOP bundle is loaded.
        let eop = IersEop::new();
        let vals = eop.eop_at(crate::J2000);
        assert_eq!(vals.dut1.value(), 0.0);
        assert_eq!(vals.xp.value(), 0.0);
    }

    #[test]
    fn iers_eop_mjd_range_is_none_without_runtime_bundle() {
        let eop = IersEop::new();
        // EOP is runtime-loaded; without a bundle the range is undefined.
        assert!(eop.mjd_range().is_none());
        assert!(eop.is_empty());
    }
}
