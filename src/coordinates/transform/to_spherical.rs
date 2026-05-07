// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cartesian-to-Spherical coordinate conversion.
//!
//! ## Scientific scope
//!
//! Spherical coordinates (longitude/latitude/distance) are the standard output
//! form for most astronomical results: star catalogues publish positions in
//! right-ascension/declination and planetary positions in ecliptic
//! longitude/latitude. Cartesian arithmetic is simpler for rotations and
//! vector operations, so the round-trip Cartesian ↔ Spherical is ubiquitous.
//!
//! ## Technical scope
//!
//! The conversion is implemented in the `affn` geometry kernel and exposed
//! through siderust's type aliases. No siderust-specific code is needed here;
//! this module exists to host the corresponding unit tests and to serve as a
//! documentation entry point.
//!
//! - `cartesian::Position::to_spherical()` uses `atan2(y, x)` for longitude
//!   and `asin(z / r)` for latitude.
//! - `spherical::Position::from_cartesian(cart)` is the complementary form.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*,
//!   4th ed. §2.1.

#[cfg(test)]
mod tests {
    use crate::coordinates::{cartesian, spherical};
    use crate::qtty::{AstronomicalUnit, Degrees};

    #[test]
    fn test_cartesian_to_spherical() {
        use crate::macros::assert_spherical_eq;
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::new(1.0, 1.0, 1.0);
        // Use affn's inherent method
        let sph: spherical::position::GCRS<AstronomicalUnit> =
            spherical::Position::from_cartesian(&cart);
        let expected = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        assert_spherical_eq!(
            &sph,
            &expected,
            1e-6,
            "Spherical coordinates do not match expected values"
        );
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        use crate::macros::assert_cartesian_eq;
        let cart_original = cartesian::position::GCRS::<AstronomicalUnit>::new(2.0, 3.0, 4.0);
        // Use affn's inherent methods
        let sph = cart_original.to_spherical();
        let cart_converted = sph.to_cartesian();
        assert_cartesian_eq!(
            &cart_original,
            &cart_converted,
            1e-6,
            "Cartesian coordinates do not match expected values"
        );
    }
}
