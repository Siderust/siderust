use super::*;
use crate::units::*;
use crate::bodies::solar_system::Moon;
use crate::astro::JulianDate;
use crate::targets::Target;
use crate::bodies::solar_system::*;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::Heliocentric, frames::Ecliptic
};

include!(concat!(env!("OUT_DIR"), "/vsop87a.rs"));

macro_rules! impl_vsop87a {
    (
        $Planet:ident,
        x: [$($x:ident),+ $(,)?],
        y: [$($y:ident),+ $(,)?],
        z: [$($z:ident),+ $(,)?]
    ) => {
        impl $Planet {
            pub fn vsop87a(jd: JulianDate) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>> {
                let (x, y, z) = position(
                    jd,
                    &[$( &$x ),+],
                    &[$( &$y ),+],
                    &[$( &$z ),+]
                );
                Target::new_static(
                    Position::new(
                        AstronomicalUnits::new(x),
                        AstronomicalUnits::new(y),
                        AstronomicalUnits::new(z)),
                    jd,
                )
            }

            pub fn vsop87a_vel(jd: JulianDate) -> Velocity<Ecliptic, AuPerDay> {
                let (vx, vy, vz) = velocity(
                    jd,
                    &[$( &$x ),+],
                    &[$( &$y ),+],
                    &[$( &$z ),+]
                );
                Velocity::new(
                    AusPerDay::new(vx),
                    AusPerDay::new(vy),
                    AusPerDay::new(vz)
                )
            }

            pub fn vsop87a_pos_vel(jd: JulianDate)
                -> (Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>>, Velocity<Ecliptic, AuPerDay>) {
                let ((x, y, z), (vx, vy, vz)) = position_velocity(
                    jd,
                    &[$( &$x ),+],
                    &[$( &$y ),+],
                    &[$( &$z ),+]
                );
                (
                    Target::new_static(Position::new(
                        AstronomicalUnits::new(x),
                        AstronomicalUnits::new(y),
                        AstronomicalUnits::new(z)), jd,),
                    Velocity::new(
                        AusPerDay::new(vx),
                        AusPerDay::new(vy),
                        AusPerDay::new(vz))
                )
            }
        }
    };
}

impl_vsop87a!(
    Mercury,
    x: [MERCURY_X0, MERCURY_X1, MERCURY_X2, MERCURY_X3, MERCURY_X4, MERCURY_X5],
    y: [MERCURY_Y0, MERCURY_Y1, MERCURY_Y2, MERCURY_Y3, MERCURY_Y4, MERCURY_Y5],
    z: [MERCURY_Z0, MERCURY_Z1, MERCURY_Z2, MERCURY_Z3, MERCURY_Z4, MERCURY_Z5]
);

impl_vsop87a!(
    Venus,
    x: [VENUS_X0, VENUS_X1, VENUS_X2, VENUS_X3, VENUS_X4, VENUS_X5],
    y: [VENUS_Y0, VENUS_Y1, VENUS_Y2, VENUS_Y3, VENUS_Y4, VENUS_Y5],
    z: [VENUS_Z0, VENUS_Z1, VENUS_Z2, VENUS_Z3, VENUS_Z4, VENUS_Z5]
);

impl_vsop87a!(
    Earth,
    x: [EARTH_X0, EARTH_X1, EARTH_X2, EARTH_X3, EARTH_X4, EARTH_X5],
    y: [EARTH_Y0, EARTH_Y1, EARTH_Y2, EARTH_Y3, EARTH_Y4, EARTH_Y5],
    z: [EARTH_Z0, EARTH_Z1, EARTH_Z2, EARTH_Z3, EARTH_Z4, EARTH_Z5]
);

impl_vsop87a!(
    Mars,
    x: [MARS_X0, MARS_X1, MARS_X2, MARS_X3, MARS_X4, MARS_X5],
    y: [MARS_Y0, MARS_Y1, MARS_Y2, MARS_Y3, MARS_Y4, MARS_Y5],
    z: [MARS_Z0, MARS_Z1, MARS_Z2, MARS_Z3, MARS_Z4, MARS_Z5]
);

impl_vsop87a!(
    Jupiter,
    x: [JUPITER_X0, JUPITER_X1, JUPITER_X2, JUPITER_X3, JUPITER_X4, JUPITER_X5],
    y: [JUPITER_Y0, JUPITER_Y1, JUPITER_Y2, JUPITER_Y3, JUPITER_Y4, JUPITER_Y5],
    z: [JUPITER_Z0, JUPITER_Z1, JUPITER_Z2, JUPITER_Z3, JUPITER_Z4, JUPITER_Z5]
);

impl_vsop87a!(
    Saturn,
    x: [SATURN_X0, SATURN_X1, SATURN_X2, SATURN_X3, SATURN_X4, SATURN_X5],
    y: [SATURN_Y0, SATURN_Y1, SATURN_Y2, SATURN_Y3, SATURN_Y4, SATURN_Y5],
    z: [SATURN_Z0, SATURN_Z1, SATURN_Z2, SATURN_Z3, SATURN_Z4, SATURN_Z5]
);

impl_vsop87a!(
    Uranus,
    x: [URANUS_X0, URANUS_X1, URANUS_X2, URANUS_X3, URANUS_X4],
    y: [URANUS_Y0, URANUS_Y1, URANUS_Y2, URANUS_Y3, URANUS_Y4],
    z: [URANUS_Z0, URANUS_Z1, URANUS_Z2, URANUS_Z3]
);

impl_vsop87a!(
    Neptune,
    x: [NEPTUNE_X0, NEPTUNE_X1, NEPTUNE_X2, NEPTUNE_X3, NEPTUNE_X4],
    y: [NEPTUNE_Y0, NEPTUNE_Y1, NEPTUNE_Y2, NEPTUNE_Y3, NEPTUNE_Y4],
    z: [NEPTUNE_Z0, NEPTUNE_Z1, NEPTUNE_Z2, NEPTUNE_Z3]
);

impl_vsop87a!(
    Moon,
    x: [EMB_X0, EMB_X1, EMB_X2, EMB_X3, EMB_X4, EMB_X5],
    y: [EMB_Y0, EMB_Y1, EMB_Y2, EMB_Y3, EMB_Y4, EMB_Y5],
    z: [EMB_Z0, EMB_Z1, EMB_Z2, EMB_Z3, EMB_Z4, EMB_Z5]
);

#[cfg(test)]
mod tests {
    use crate::units::AU;
    use crate::astro::JulianDate;
    use super::*;
    use crate::coordinates::cartesian::Position;
    use crate::macros::assert_cartesian_eq;

    const PRECISION: f64 = 1e-3;

    /// Test Mercury's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_mercury_at_epoch() {

        let coord = Mercury::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(-0.1302524*AU, -0.4472397*AU, -0.0245799*AU), PRECISION);
    }

    /// Test Venus heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_venus_at_epoch() {
        // At epoch, compute heliocentric coordinates
        let coord = Venus::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(-0.7183022991131299*AU, -0.03265428553900499*AU, 0.040809*AU), PRECISION);
    }

    /// Test Mars's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_mars_at_epoch() {
        // At epoch, compute heliocentric coordinates
        let coord = Mars::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(1.3907159447538169*AU, -0.013416322699311728*AU, -0.034668*AU), PRECISION);
    }

    /// Test Jupiter's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_jupiter_at_epoch() {
        // At epoch, compute heliocentric coordinates
        let coord = Jupiter::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(4.008895*AU, 2.940636*AU, -0.101869*AU), 1e-2);
    }

    /// Test Saturn's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_saturn_at_epoch() {
        // At epoch, compute heliocentric coordinates
        let coord = Saturn::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(6.412182*AU, 6.572783*AU, -0.369816*AU), 5e-2);
    }

    /// Test Uranus's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_uranus_at_epoch() {
        // At epoch, compute heliocentric coordinates
        let coord = Uranus::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(14.438269*AU, -13.733294*AU, -0.238515*AU), 1e-2);
    }

    /// Test Neptune's heliocentric coordinates at epoch J2000.0
    #[test]
    fn test_neptune_at_epoch() {
        let coord = Neptune::vsop87a(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(coord, Position::new(16.817474*AU, -24.990018*AU, 0.126993*AU), 1e-2);
    }
}
