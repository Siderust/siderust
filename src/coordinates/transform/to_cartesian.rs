
// Note: Spherical-to-Cartesian conversion is handled by affn's inherent methods:
// - spherical::Position::to_cartesian()
// - cartesian::Position::from_spherical()
//
// These are re-exported through siderust's type aliases, so users can call:
//   sph_position.to_cartesian()  
//   cartesian::Position::from_spherical(&sph_position)

#[cfg(test)]
mod tests {
    use crate::coordinates::{cartesian, spherical};
    use qtty::{AstronomicalUnit, Degrees};

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
