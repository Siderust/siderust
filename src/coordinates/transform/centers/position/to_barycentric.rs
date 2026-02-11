// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::JulianDate;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::CenterShiftProvider;
use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{self, MutableFrame},
};
use qtty::{LengthUnit, Quantity};

// =============================================================================
// Heliocentric → Barycentric (pure translation)
// =============================================================================
//
// This implementation delegates to CenterShiftProvider, which sources positions
// from the ephemeris trait. Using default AstroContext provides VSOP87 positions.

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<qtty::AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Heliocentric, Barycentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Heliocentric, Barycentric, F>>::shift(jd, &ctx);

        // Apply the shift
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Barycentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
    }
}

// =============================================================================
// Geocentric → Barycentric (pure translation)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<qtty::AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Geocentric, Barycentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Geocentric, Barycentric, F>>::shift(jd, &ctx);

        // Apply the shift
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Barycentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
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
