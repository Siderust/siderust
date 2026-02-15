// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Orientation Parameters (EOP)
//!
//! Earth Orientation Parameters connect the celestial (GCRS) and terrestrial
//! (ITRS) reference frames by accounting for irregularities in Earth's rotation.
//!
//! The key EOP quantities are:
//!
//! | Parameter | Description |
//! |-----------|-------------|
//! | **UT1 − UTC** | Difference between Universal Time and Coordinated Universal Time |
//! | **xₚ, yₚ** | Pole coordinates — position of the CIP in the ITRS |
//! | **dX, dY** | Celestial pole offsets — corrections to the IAU precession-nutation model |
//! | **LOD** | Length Of Day excess (deviation from 86400 SI seconds) |
//!
//! ## Architecture
//!
//! The [`EopProvider`] trait abstracts over different EOP data sources:
//! - [`NullEop`]: zero-correction placeholder (assumes UTC ≈ UT1, no polar motion)
//! - Future: IERS Bulletin A/B parsers, finals2000A, etc.
//!
//! ## References
//!
//! * IERS Conventions (2010), Chapter 5
//! * IERS Technical Note 36

use crate::time::JulianDate;
use qtty::*;

/// A set of Earth Orientation Parameters at a given epoch.
///
/// Field units match IERS publications for natural construction from
/// Bulletin A/B data:
/// - `dut1` in [`Seconds`] (UT1 − UTC)
/// - `xp`, `yp` in [`Arcseconds`] (pole coordinates)
/// - `dx`, `dy` in [`MilliArcseconds`] (celestial pole offsets)
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
}

impl Default for EopValues {
    fn default() -> Self {
        Self {
            dut1: Seconds::new(0.0),
            xp: Arcseconds::new(0.0),
            yp: Arcseconds::new(0.0),
            dx: MilliArcseconds::new(0.0),
            dy: MilliArcseconds::new(0.0),
        }
    }
}

impl EopValues {
    /// Convert UTC Julian Date to UT1 Julian Date using this EOP's dUT1.
    #[inline]
    pub fn jd_ut1(&self, jd_utc: JulianDate) -> JulianDate {
        JulianDate::new(jd_utc.value() + self.dut1.to::<Day>().value())
    }
}

/// Trait for providing Earth Orientation Parameters at a given epoch.
///
/// Implementors might read from IERS data files, interpolate tables,
/// or return zero corrections for simplified calculations.
pub trait EopProvider {
    /// Look up (or interpolate) EOP values for the given Julian Date (UTC).
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
    fn eop_at(&self, _jd_utc: JulianDate) -> EopValues {
        EopValues::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_eop_returns_zeros() {
        let eop = NullEop;
        let vals = eop.eop_at(JulianDate::J2000);
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
        let jd_utc = JulianDate::J2000;
        let jd_ut1 = vals.jd_ut1(jd_utc);
        let diff_days = (jd_ut1 - jd_utc).value();
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
}
