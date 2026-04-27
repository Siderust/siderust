// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomical time units not already present in `qtty`.
//!
//! [`qtty`] ships [`crate::qtty::time::SiderealYear`] (365.256 363 004 d, the actual
//! measured period of Earth relative to the fixed stars) and
//! [`crate::qtty::time::JulianYear`] (365.25 d exactly).
//!
//! This module adds the **Gaussian year** — the year length that is implied by
//! the Gaussian gravitational constant `k = 0.01720209895 AU^{3/2} d^{-1}`.
//! For a test particle at 1 AU from the Sun, Kepler's third law gives:
//!
//! ```text
//! T_Gaussian = 2π / k ≈ 365.256 898 326 d
//! ```
//!
//! This is the natural "year" unit for expressing orbital periods in the
//! heliocentric AU-day system: `T (Gaussian years) = a^{3/2} (AU)`.
//! It differs from the measured sidereal year by ≈ 46 s, because the
//! Gaussian constant encodes the solar mass + test particle mass, while
//! the modern sidereal year accounts for the full two-body problem with
//! a more precise solar GM.

use crate::qtty::{Quantity, Time, Unit};

// ─────────────────────────────────────────────────────────────────────────────
// GaussianYear — unit definition
// ─────────────────────────────────────────────────────────────────────────────

/// The **Gaussian year**: the orbital period implied by the Gaussian
/// gravitational constant `k = 0.01720209895 AU^{3/2} d^{-1}`.
///
/// ```text
/// T_Gaussian = 2π / k ≈ 365.256 898 326 d  ≈  31 558 196.0 s
/// ```
///
/// This is the natural time unit for Kepler's third law in the AU-day system:
///
/// ```text
/// T [Gaussian years]  =  a [AU]^{3/2}
/// ```
///
/// Use [`GaussianYears`] for values and `.to::<Day>()` (or any other
/// [`qtty`] time unit) for unit-safe conversion.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct GaussianYear;

/// Seconds in one Gaussian year (`2π / k * 86 400`, where `k = 0.01720209895`).
pub(crate) const GAUSSIAN_YEAR_SECONDS: f64 = 365.256_898_326 * 86_400.0;

impl Unit for GaussianYear {
    /// Conversion factor to the canonical time unit (seconds).
    const RATIO: f64 = GAUSSIAN_YEAR_SECONDS;
    type Dim = Time;
    const SYMBOL: &'static str = "gyr";
}

/// A quantity measured in Gaussian years.
///
/// Supports `.to::<Day>()`, `.to::<Second>()`, etc., via the standard `qtty`
/// unit-conversion machinery.
///
/// # Example
///
/// ```rust
/// use siderust::astro::units::GaussianYears;
/// use siderust::qtty::time::Day;
///
/// // Earth's orbital period via Kepler's 3rd law (a = 1 AU → T = 1 Gaussian year).
/// let period_days = GaussianYears::new(1.0).to::<Day>();
/// assert!((period_days.value() - 365.256_898_326).abs() < 1e-9);
/// ```
pub type GaussianYears = Quantity<GaussianYear>;

/// One Gaussian year as a typed constant.
pub const GAUSSIAN_YEAR: GaussianYears = GaussianYears::new(1.0);

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::time::{Day, JulianYear, Second, SiderealYear};

    #[test]
    fn gaussian_year_to_days() {
        let d = GaussianYears::new(1.0).to::<Day>();
        assert!((d.value() - 365.256_898_326).abs() < 1e-9);
    }

    #[test]
    fn gaussian_year_to_seconds() {
        let s = GaussianYears::new(1.0).to::<Second>();
        assert!((s.value() - GAUSSIAN_YEAR_SECONDS).abs() < 1e-3);
    }

    #[test]
    fn gaussian_year_differs_from_sidereal_year() {
        // Should differ by ~46 s (≈ 0.000535 d).
        let g = GaussianYears::new(1.0).to::<Day>().value();
        let s = crate::qtty::Quantity::<SiderealYear>::new(1.0).to::<Day>().value();
        let diff = (g - s).abs();
        assert!(diff > 4e-4 && diff < 7e-4, "diff = {diff}");
    }

    #[test]
    fn gaussian_year_differs_from_julian_year() {
        // Julian year = 365.25 d; Gaussian year ≈ 365.257 d.
        let g = GaussianYears::new(1.0).to::<Day>().value();
        let j = crate::qtty::Quantity::<JulianYear>::new(1.0).to::<Day>().value();
        assert!((g - j).abs() > 5e-3);
    }

    #[test]
    fn roundtrip_gaussian_year_seconds() {
        let original = GaussianYears::new(3.7);
        let seconds = original.to::<Second>();
        let back = seconds.to::<GaussianYear>();
        assert!((back.value() - original.value()).abs() < 1e-12);
    }

    #[test]
    fn gaussian_year_symbol() {
        assert_eq!(GaussianYear::SYMBOL, "gyr");
    }
}
