
// Note: Cartesian-to-Spherical conversion is handled by affn's inherent methods:
// - cartesian::Position::to_spherical()
// - spherical::Position::from_cartesian()
//
// These are re-exported through siderust's type aliases, so users can call:
//   cart_position.to_spherical()
//   spherical::Position::from_cartesian(&cart_position)

#[cfg(test)]
mod tests {
    use crate::coordinates::spherical::ext::IcrsPositionExt;
    use crate::coordinates::{cartesian, spherical};
    use qtty::{AstronomicalUnit, Degrees};

    #[test]
    fn test_cartesian_to_spherical() {
        use crate::macros::assert_spherical_eq;
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::new(1.0, 1.0, 1.0);
        // Use affn's inherent method
        let sph: spherical::position::GCRS<AstronomicalUnit> = cart.to_spherical();
        let expected = spherical::position::GCRS::<AstronomicalUnit>::new_icrs(
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
