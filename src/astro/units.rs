// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical Units (siderust extensions)
//!
//! Astronomical time units that are not already shipped by [`qtty`], in
//! particular the **Gaussian year** used to express orbital periods in the
//! heliocentric AU-day system.
//!
//! ## Scientific scope
//!
//! Two "year" units are widely used in astronomy: the **sidereal year**
//! (≈ 365.256 363 004 d, the actual measured period of Earth relative to
//! the fixed stars) and the **Julian year** (365.25 d exactly, used for
//! Julian centuries). Both are provided by `qtty`. This module adds the
//! **Gaussian year**, the period implied by Gauss's gravitational constant
//! `k = 0.01720209895 AU^{3/2} d^{-1}`. It is the natural time unit for
//! Kepler's third law in the AU-day system, since for a test particle at
//! 1 AU `T_Gaussian = 2π / k ≈ 365.256 898 326 d` and
//! `T [gyr] = a [AU]^{3/2}`. It differs from the sidereal year by ≈ 46 s
//! because `k` encodes only the Sun + test-particle two-body problem.
//!
//! ## Technical scope
//!
//! The unit is implemented as a zero-sized `GaussianYear` type implementing
//! the [`Unit`] trait of [`qtty`], with `RATIO = 365.256 898 326 × 86 400`
//! seconds and dimension [`Time`]. The corresponding quantity alias
//! [`GaussianYears`] supports the standard `qtty` conversion machinery
//! (`.to::<Day>()`, `.to::<Second>()`, …) and is interoperable with every
//! other time unit in the crate.
//!
//! ## References
//!
//! * Gauss, C. F., *Theoria Motus Corporum Coelestium* (1809)
//! * Seidelmann, *Explanatory Supplement to the Astronomical Almanac*,
//!   §8 (time units)
//! * IAU 2012 Resolution B2 (redefinition of the astronomical unit)

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
// Heliocentric period helper
// ─────────────────────────────────────────────────────────────────────────────

/// Heliocentric orbital period in days from Kepler's third law (AU-day system).
///
/// Applies `T [Gaussian years] = a [AU]^{3/2}` and converts to days using the
/// Gaussian year (2π / k ≈ 365.256 898 326 d). Valid only for heliocentric
/// (Sun-centered) orbits; for other central bodies derive mean motion from the
/// body's gravitational parameter instead.
///
/// # Arguments
///
/// - `semi_major_axis_au` — semi-major axis in astronomical units.
#[inline]
pub(crate) fn heliocentric_period_days(semi_major_axis_au: f64) -> f64 {
    use crate::qtty::time::Day;
    GaussianYears::new(semi_major_axis_au * semi_major_axis_au.sqrt())
        .to::<Day>()
        .value()
}

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
        let s = crate::qtty::Quantity::<SiderealYear>::new(1.0)
            .to::<Day>()
            .value();
        let diff = (g - s).abs();
        assert!(diff > 4e-4 && diff < 7e-4, "diff = {diff}");
    }

    #[test]
    fn gaussian_year_differs_from_julian_year() {
        // Julian year = 365.25 d; Gaussian year ≈ 365.257 d.
        let g = GaussianYears::new(1.0).to::<Day>().value();
        let j = crate::qtty::Quantity::<JulianYear>::new(1.0)
            .to::<Day>()
            .value();
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
