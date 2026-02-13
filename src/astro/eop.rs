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
#[derive(Debug, Clone, Copy)]
pub struct EopValues {
    /// UT1 − UTC in seconds.
    pub dut1: f64,
    /// Pole x-coordinate (arcseconds).
    pub xp: f64,
    /// Pole y-coordinate (arcseconds).
    pub yp: f64,
    /// Celestial pole offset dX (milliarcseconds).
    pub dx: f64,
    /// Celestial pole offset dY (milliarcseconds).
    pub dy: f64,
}

impl Default for EopValues {
    fn default() -> Self {
        Self {
            dut1: 0.0,
            xp: 0.0,
            yp: 0.0,
            dx: 0.0,
            dy: 0.0,
        }
    }
}

impl EopValues {
    /// Convert UTC Julian Date to UT1 Julian Date using this EOP's dUT1.
    #[inline]
    pub fn jd_ut1(&self, jd_utc: JulianDate) -> JulianDate {
        JulianDate::new(jd_utc.value() + self.dut1 / 86400.0)
    }

    /// Pole x-coordinate in radians.
    #[inline]
    pub fn xp_rad(&self) -> f64 {
        self.xp * std::f64::consts::PI / (180.0 * 3600.0)
    }

    /// Pole y-coordinate in radians.
    #[inline]
    pub fn yp_rad(&self) -> f64 {
        self.yp * std::f64::consts::PI / (180.0 * 3600.0)
    }

    /// Celestial pole offset dX in radians.
    #[inline]
    pub fn dx_rad(&self) -> f64 {
        self.dx * 1e-3 * std::f64::consts::PI / (180.0 * 3600.0)
    }

    /// Celestial pole offset dY in radians.
    #[inline]
    pub fn dy_rad(&self) -> f64 {
        self.dy * 1e-3 * std::f64::consts::PI / (180.0 * 3600.0)
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
        assert_eq!(vals.dut1, 0.0);
        assert_eq!(vals.xp, 0.0);
        assert_eq!(vals.yp, 0.0);
        assert_eq!(vals.dx, 0.0);
        assert_eq!(vals.dy, 0.0);
    }

    #[test]
    fn eop_jd_ut1_conversion() {
        let vals = EopValues {
            dut1: 0.35, // 350 ms ahead
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
            xp: 0.1,  // 0.1″
            yp: 0.2,  // 0.2″
            dx: 0.3,  // 0.3 mas
            dy: 0.4,  // 0.4 mas
            ..Default::default()
        };
        assert!((vals.xp_rad() - 0.1 * as2rad).abs() < 1e-20);
        assert!((vals.yp_rad() - 0.2 * as2rad).abs() < 1e-20);
        assert!((vals.dx_rad() - 0.3e-3 * as2rad).abs() < 1e-20);
        assert!((vals.dy_rad() - 0.4e-3 * as2rad).abs() < 1e-20);
    }
}
