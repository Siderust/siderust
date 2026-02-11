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
// Geocentric → Heliocentric (pure translation, no aberration)
// =============================================================================
//
// Center transforms are pure geometry: they only translate positions from one
// origin to another. Aberration is an observation-model effect that depends on
// observer velocity and must be applied explicitly via the observation module.
//
// This implementation delegates to CenterShiftProvider, which sources positions
// from the ephemeris trait. Using default AstroContext provides VSOP87 positions.

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<qtty::AstronomicalUnits> + PartialEq + std::fmt::Debug,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Geocentric, Heliocentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Geocentric, Heliocentric, F>>::shift(jd, &ctx);

        // Apply the shift
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Heliocentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
    }
}

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<qtty::AstronomicalUnits> + PartialEq + std::fmt::Debug,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Barycentric, Heliocentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Barycentric, Heliocentric, F>>::shift(jd, &ctx);

        // Apply the shift
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Heliocentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Sun;
    use crate::coordinates::cartesian::position::Ecliptic;
    use crate::coordinates::centers::*;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use qtty::Au;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = *Sun::vsop87e(JulianDate::J2000).get_position();
        let sun_helio: Ecliptic<Au, Heliocentric> = sun_bary.transform(JulianDate::J2000);
        let expected_sun_helio = Ecliptic::<Au>::CENTER;
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            EPSILON,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }

    #[test] // Geocentric -> Heliocentric
    fn test_geo() {
        let sun = *Sun::vsop87e(JulianDate::J2000).get_position();
        let sun_geo: Ecliptic<Au, Geocentric> = sun.transform(JulianDate::J2000);
        let sun_helio: Ecliptic<Au> = sun_geo.transform(JulianDate::J2000);
        let expected_sun_helio = Ecliptic::<Au>::CENTER;
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            1e-8,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }
}
