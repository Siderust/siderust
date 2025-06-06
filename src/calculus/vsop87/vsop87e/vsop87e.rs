use crate::bodies::solar_system::*;
use crate::calculus::vsop87::compute_vsop87;
use crate::coordinates::{cartesian::Position, centers::Barycentric, frames::Ecliptic};
use crate::targets::Target;
use crate::units::JulianDay;
use super::{sun_vsop87e::*,
            mercury_vsop87e::*,
            venus_vsop87e::*,
            earth_vsop87e::*,
            mars_vsop87e::*,
            jupiter_vsop87e::*,
            saturn_vsop87e::*,
            uranus_vsop87e::*,
            neptune_vsop87e::*};

macro_rules! impl_vsop87e {
    (
        $Planet:ident,
        x: [$($x:ident),+ $(,)?],
        y: [$($y:ident),+ $(,)?],
        z: [$($z:ident),+ $(,)?]
    ) => {
        impl $Planet {
            pub fn vsop87e(jd: JulianDay) -> Target<Position<Barycentric, Ecliptic>> {
                let (x, y, z) = compute_vsop87(
                    jd,
                    &[$( &$x ),+],
                    &[$( &$y ),+],
                    &[$( &$z ),+]
                );
                Target::new_static(
                    Position::<Barycentric, Ecliptic>::new(x, y, z),
                    jd,
                )
            }
        }
    };
}

impl Sun {
    pub fn vsop87e(jd: JulianDay) -> Target<Position<Barycentric, Ecliptic>> {
        let (x, y, z) = compute_vsop87(
            jd,
            &[&SUN_X0, &SUN_X1, &SUN_X2, &SUN_X3, &SUN_X4, &SUN_X5],
            &[&SUN_Y0, &SUN_Y1, &SUN_Y2, &SUN_Y3, &SUN_Y4, &SUN_Y5],
            &[&SUN_Z0, &SUN_Z1, &SUN_Z2, &SUN_Z3, &SUN_Z4, &SUN_Z5]
        );
        Target::new_static(
            Position::<Barycentric, Ecliptic>::new(x, y, z),
            jd,
        )
    }
}

impl_vsop87e!(
    Mercury,
    x: [MERCURY_X0, MERCURY_X1, MERCURY_X2, MERCURY_X3, MERCURY_X4, MERCURY_X5],
    y: [MERCURY_Y0, MERCURY_Y1, MERCURY_Y2, MERCURY_Y3, MERCURY_Y4, MERCURY_Y5],
    z: [MERCURY_Z0, MERCURY_Z1, MERCURY_Z2, MERCURY_Z3, MERCURY_Z4, MERCURY_Z5]
);

impl_vsop87e!(
    Venus,
    x: [VENUS_X0, VENUS_X1, VENUS_X2, VENUS_X3, VENUS_X4, VENUS_X5],
    y: [VENUS_Y0, VENUS_Y1, VENUS_Y2, VENUS_Y3, VENUS_Y4, VENUS_Y5],
    z: [VENUS_Z0, VENUS_Z1, VENUS_Z2, VENUS_Z3, VENUS_Z4, VENUS_Z5]
);

impl_vsop87e!(
    Earth,
    x: [EARTH_X0, EARTH_X1, EARTH_X2, EARTH_X3, EARTH_X4, EARTH_X5],
    y: [EARTH_Y0, EARTH_Y1, EARTH_Y2, EARTH_Y3, EARTH_Y4, EARTH_Y5],
    z: [EARTH_Z0, EARTH_Z1, EARTH_Z2, EARTH_Z3, EARTH_Z4, EARTH_Z5]
);

impl_vsop87e!(
    Mars,
    x: [MARS_X0, MARS_X1, MARS_X2, MARS_X3, MARS_X4, MARS_X5],
    y: [MARS_Y0, MARS_Y1, MARS_Y2, MARS_Y3, MARS_Y4, MARS_Y5],
    z: [MARS_Z0, MARS_Z1, MARS_Z2, MARS_Z3, MARS_Z4, MARS_Z5]
);

impl_vsop87e!(
    Jupiter,
    x: [JUPITER_X0, JUPITER_X1, JUPITER_X2, JUPITER_X3, JUPITER_X4, JUPITER_X5],
    y: [JUPITER_Y0, JUPITER_Y1, JUPITER_Y2, JUPITER_Y3, JUPITER_Y4, JUPITER_Y5],
    z: [JUPITER_Z0, JUPITER_Z1, JUPITER_Z2, JUPITER_Z3, JUPITER_Z4, JUPITER_Z5]
);

impl_vsop87e!(
    Saturn,
    x: [SATURN_X0, SATURN_X1, SATURN_X2, SATURN_X3, SATURN_X4, SATURN_X5],
    y: [SATURN_Y0, SATURN_Y1, SATURN_Y2, SATURN_Y3, SATURN_Y4, SATURN_Y5],
    z: [SATURN_Z0, SATURN_Z1, SATURN_Z2, SATURN_Z3, SATURN_Z4, SATURN_Z5]
);

impl_vsop87e!(
    Uranus,
    x: [URANUS_X0, URANUS_X1, URANUS_X2, URANUS_X3, URANUS_X4],
    y: [URANUS_Y0, URANUS_Y1, URANUS_Y2, URANUS_Y3, URANUS_Y4],
    z: [URANUS_Z0, URANUS_Z1, URANUS_Z2, URANUS_Z3]
);

impl_vsop87e!(
    Neptune,
    x: [NEPTUNE_X0, NEPTUNE_X1, NEPTUNE_X2, NEPTUNE_X3, NEPTUNE_X4],
    y: [NEPTUNE_Y0, NEPTUNE_Y1, NEPTUNE_Y2, NEPTUNE_Y3, NEPTUNE_Y4],
    z: [NEPTUNE_Z0, NEPTUNE_Z1, NEPTUNE_Z2, NEPTUNE_Z3]
);

/*
#[cfg(test)]
mod tests {
    use crate::{time::julian_date::J2000};
    use crate::coordinates::{Position, frames::Ecliptic, centers::Heliocentric};

    /// Helper function to compare two floating-point numbers with a tolerance.
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn check_cartesian(actual: Position<Heliocentric, Ecliptic>, expected_x: f64, expected_y: f64, expected_z: f64, tol: f64) {
        assert!(approx_eq(actual.x(), expected_x, tol), "current x = {}, expected x = {}", actual.x(), expected_x);
        assert!(approx_eq(actual.y(), expected_y, tol), "current y = {}, expected y = {}", actual.y(), expected_y);
        assert!(approx_eq(actual.z(), expected_z, tol), "current z = {}, expected z = {}", actual.z(), expected_z);
    }

    /// Test Mercury's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_mercury_at_epoch() {
        use crate::bodies::Mercury;

        let coord = Mercury::vsop87e(J2000);
        check_cartesian(coord, -0.1302524, -0.4472397, -0.0245799, 1e-3);
    }

    /// Test Venus heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_venus_at_epoch() {
        use crate::bodies::Venus;

        // At epoch, compute heliocentric coordinates
        let coord = Venus::vsop87e(J2000);
        check_cartesian(coord, -0.7183022991131299, -0.03265428553900499, 0.040809, 1e-3);
    }

    /// Test Mars's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_mars_at_epoch() {
        use crate::bodies::Mars;

        // At epoch, compute heliocentric coordinates
        let coord = Mars::vsop87e(J2000);
        check_cartesian(coord, 1.3907159447538169, -0.013416322699311728, -0.034668, 1e-3);
    }

    /// Test Jupiter's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_jupiter_at_epoch() {
        use crate::bodies::Jupiter;

        // At epoch, compute heliocentric coordinates
        let coord = Jupiter::vsop87e(J2000);
        check_cartesian(coord, 4.008895, 2.940636, -0.101869, 1e-2);
    }

    /// Test Saturn's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_saturn_at_epoch() {
        use crate::bodies::Saturn;

        // At epoch, compute heliocentric coordinates
        let coord = Saturn::vsop87e(J2000);
        check_cartesian(coord, 6.412182, 6.572783, -0.369816, 5e-2);
    }

    /// Test Uranus's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_uranus_at_epoch() {
        use crate::bodies::Uranus;

        // At epoch, compute heliocentric coordinates
        let coord = Uranus::vsop87e(J2000);
        check_cartesian(coord, 14.438269, -13.733294, -0.238515, 1e-2);
    }

    /// Test Neptune's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_neptune_at_epoch() {
        use crate::bodies::Neptune;

        let coord = Neptune::vsop87e(J2000);
        check_cartesian(coord, 16.817474, -24.990018, 0.126993, 1e-2);
    }
}
*/