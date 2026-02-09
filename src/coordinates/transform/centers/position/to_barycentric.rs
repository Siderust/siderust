// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::{Earth, Sun};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::position::Ecliptic,
    cartesian::Position,
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{self, MutableFrame},
    transform::Transform,
};
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// Heliocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let sun_bary_ecl_au = *Sun::vsop87e(jd).get_position();

        // VSOP87 gives the Sun's position in AstronomicalUnits
        let sun_bary_ecl = Ecliptic::<U, Barycentric>::new(
            sun_bary_ecl_au.x(),
            sun_bary_ecl_au.y(),
            sun_bary_ecl_au.z(),
        );

        let sun_bary_f: Position<Barycentric, F, U> = sun_bary_ecl.transform(jd); // (Bary-Ecl) -> (Bary-F)
        Position::from_vec3_origin(self.as_vec3() + sun_bary_f.as_vec3())
    }
}

// Geocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let earth_bary_ecl_au = *Earth::vsop87e(jd).get_position();

        // VSOP87 gives the Earth's position in AstronomicalUnits
        let earth_bary_ecl = Ecliptic::<U, Barycentric>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        let earth_bary_f: Position<Barycentric, F, U> = earth_bary_ecl.transform(jd);
        Position::from_vec3_origin(self.as_vec3() + earth_bary_f.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::{Earth, Sun};
    use crate::coordinates::cartesian::position::Ecliptic;
    use crate::coordinates::centers::*;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use qtty::Au;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Heliocentric -> Barycentric
    fn test_helio() {
        let sun_helio = Ecliptic::<Au>::CENTER;
        let sun_bary: Ecliptic<Au, Barycentric> = sun_helio.transform(JulianDate::J2000);
        let expected_sun_bary = *Sun::vsop87e(JulianDate::J2000).get_position();
        assert_cartesian_eq!(sun_bary, expected_sun_bary, EPSILON);
    }

    #[test] // Geocentric -> Barycentric
    fn test_geo() {
        let earth_geo = Ecliptic::<Au, Geocentric>::CENTER;
        let earth_bary: Ecliptic<Au, Barycentric> = earth_geo.transform(JulianDate::J2000);
        let expected_earth_bary = *Earth::vsop87e(JulianDate::J2000).get_position();
        assert_cartesian_eq!(&earth_bary, &expected_earth_bary, EPSILON);
    }
}
