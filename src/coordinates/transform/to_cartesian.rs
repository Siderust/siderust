// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical-to-Cartesian coordinate conversion.
//!
//! ## Scientific scope
//!
//! Spherical coordinates (longitude/latitude or right-ascension/declination
//! plus distance) are natural for observational astronomy, while Cartesian
//! coordinates are natural for linear algebra (rotations, vector arithmetic,
//! distance calculations). The conversion between the two representations is
//! lossless (up to floating-point precision).
//!
//! ## Technical scope
//!
//! The conversion is implemented in the `affn` geometry kernel and exposed
//! through siderust's type aliases. No siderust-specific code is needed here;
//! this module exists to host the corresponding unit tests and to serve as a
//! documentation entry point.
//!
//! - `spherical::Position::to_cartesian()` applies:
//!   `x = r cos(lat) cos(lon)`, `y = r cos(lat) sin(lon)`, `z = r sin(lat)`.
//! - `cartesian::Position::from_spherical(sph)` is the complementary entry point.
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
    fn test_spherical_to_cartesian() {
        use crate::macros::assert_cartesian_eq;
        let sph = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        // Use affn's inherent method
        let cart = sph.to_cartesian();
        let expected = cartesian::position::GCRS::new(1.0, 1.0, 1.0);
        assert_cartesian_eq!(
            &cart,
            &expected,
            1e-6,
            "Cartesian coordinates do not match expected values"
        );
    }

    #[test]
    fn test_spherical_cartesian_round_trip() {
        use crate::macros::assert_spherical_eq;
        let sph_original = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
            5.0,
        );
        // Use affn's inherent methods
        let cart = sph_original.to_cartesian();
        let sph_converted = spherical::Position::from_cartesian(&cart);
        assert_spherical_eq!(
            &sph_original,
            &sph_converted,
            1e-6,
            "Spherical coordinates do not match expected values"
        );
    }
}
