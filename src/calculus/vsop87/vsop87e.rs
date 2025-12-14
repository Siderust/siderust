use super::*;
use crate::astro::JulianDate;
use crate::bodies::solar_system::*;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::Barycentric,
    frames::Ecliptic,
};
use crate::targets::Target;
use qtty::*;

#[allow(clippy::approx_constant)]
#[rustfmt::skip]
mod vsop_data {
    include!(concat!(env!("OUT_DIR"), "/vsop87e.rs"));
}
use vsop_data::*;

macro_rules! impl_vsop87e {
    (
        $Planet:ident,
        x: [$($x:ident),+ $(,)?],
        y: [$($y:ident),+ $(,)?],
        z: [$($z:ident),+ $(,)?]
    ) => {
        impl $Planet {
            pub fn vsop87e(jd: JulianDate) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>> {
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

            pub fn vsop87e_vel(jd: JulianDate) -> Velocity<Ecliptic, AuPerDay> {
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

            pub fn vsop87e_pos_vel(jd: JulianDate)
                -> (Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>, Velocity<Ecliptic, AuPerDay>) {
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

impl Sun {
    pub fn vsop87e(jd: JulianDate) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>> {
        let (x, y, z) = position(
            jd,
            &[&SUN_X0, &SUN_X1, &SUN_X2, &SUN_X3, &SUN_X4, &SUN_X5],
            &[&SUN_Y0, &SUN_Y1, &SUN_Y2, &SUN_Y3, &SUN_Y4, &SUN_Y5],
            &[&SUN_Z0, &SUN_Z1, &SUN_Z2, &SUN_Z3, &SUN_Z4, &SUN_Z5],
        );
        Target::new_static(
            Position::new(
                AstronomicalUnits::new(x),
                AstronomicalUnits::new(y),
                AstronomicalUnits::new(z),
            ),
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

#[cfg(test)]
mod tests {
    use crate::astro::JulianDate;
    use crate::coordinates::cartesian::Position;
    use crate::macros::assert_cartesian_eq;
    use qtty::AU;

    const PRECISION: f64 = 1.0e-6;

    /// Mercury barycentric (ecliptic J2000) at epoch J2000.0 (JD_TDB = 2451545.0)
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_mercury_at_epoch() {
        use crate::bodies::Mercury;
        let coord = Mercury::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(-0.1372349394 * AU, -0.4500758422 * AU, -0.0243922379 * AU),
            PRECISION
        );
    }

    /// Venus barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_venus_at_epoch() {
        use crate::bodies::Venus;
        let coord = Venus::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(-0.7254438061 * AU, -0.0354427729 * AU, 0.0412204390 * AU),
            PRECISION
        );
    }

    /// Test Earth's barycentric coordinates at epoch J2000.0 (VSOP87E)
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_earth_barycentric_at_epoch() {
        use crate::bodies::Earth;
        let coord = Earth::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(-0.1842769894 * AU, 0.9644534522 * AU, 0.0002021000 * AU),
            PRECISION
        );
    }

    /// Mars barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_mars_at_epoch() {
        use crate::bodies::Mars;
        let coord = Mars::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(1.3835744053 * AU, -0.0162038666 * AU, -0.0342616574 * AU),
            PRECISION
        );
    }

    /// Jupiter barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_jupiter_at_epoch() {
        use crate::bodies::Jupiter;
        let coord = Jupiter::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(3.9940325025 * AU, 2.9357928287 * AU, -0.1015776146 * AU),
            PRECISION
        );
    }

    /// Saturn barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_saturn_at_epoch() {
        use crate::bodies::Saturn;
        let coord = Saturn::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(6.3992653459 * AU, 6.5672047374 * AU, -0.3688706482 * AU),
            PRECISION
        );
    }

    /// Uranus barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_uranus_at_epoch() {
        use crate::bodies::Uranus;
        let coord = Uranus::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(14.4247519568 * AU, -13.7371045087 * AU, -0.2379360887 * AU),
            PRECISION
        );
    }

    /// Neptune barycentric (ecliptic J2000) at epoch J2000.0
    #[test]
    #[cfg_attr(siderust_stubs, ignore)]
    fn test_neptune_at_epoch() {
        use crate::bodies::Neptune;
        let coord = Neptune::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(
            coord,
            Position::new(16.8049701269 * AU, -24.9944513569 * AU, 0.1274251215 * AU),
            PRECISION
        );
    }
}
