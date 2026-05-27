// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Eclipse classification for solar/shadow geometry.

use affn::cartesian::Position;
use affn::{ReferenceCenter, ReferenceFrame};
use qtty::length::{Meter, Meters};

/// Eclipse classification for a satellite illuminated by a single light
/// source (typically the Sun) with one obscuring body (typically Earth).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EclipseState {
    /// Satellite is fully illuminated.
    Sunlight,
    /// Penumbra: light source is partially obscured.
    Penumbra,
    /// Umbra: light source is fully obscured.
    Umbra,
}

/// Classify the eclipse state of a satellite illuminated by a Sun-like
/// light source and obscured by a single occulting body.
///
/// All positions must be expressed in the same Cartesian center and frame.
/// `sat` is the satellite position, `sun` and `occulter_center` are the
/// positions of the light source and the obscuring body, and the radii
/// are the physical body radii.
pub fn eclipse_state<C, F>(
    sat: &Position<C, F, Meter>,
    sun: &Position<C, F, Meter>,
    sun_radius: Meters,
    occulter_center: &Position<C, F, Meter>,
    occulter_radius: Meters,
) -> EclipseState
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    // Reduce to "Sun seen from satellite vs occulter seen from satellite".
    let f = super::occultation::occultation_fraction(
        sat,
        sun,
        sun_radius,
        occulter_center,
        occulter_radius,
    );
    if f >= 1.0 - 1e-9 {
        EclipseState::Umbra
    } else if f > 0.0 {
        EclipseState::Penumbra
    } else {
        EclipseState::Sunlight
    }
}

/// True iff the satellite is in any kind of eclipse (penumbra or umbra).
pub fn solar_eclipsing<C, F>(
    sat: &Position<C, F, Meter>,
    sun: &Position<C, F, Meter>,
    sun_radius: Meters,
    occulter_center: &Position<C, F, Meter>,
    occulter_radius: Meters,
) -> bool
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    !matches!(
        eclipse_state(sat, sun, sun_radius, occulter_center, occulter_radius),
        EclipseState::Sunlight
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use affn::cartesian::Position;
    use affn::DeriveReferenceCenter;
    use affn::DeriveReferenceFrame;
    use qtty::length::{Meter, Meters};
    use qtty::Quantity;

    #[derive(Debug, Clone, Copy, DeriveReferenceCenter)]
    struct C;

    #[derive(Debug, Clone, Copy, DeriveReferenceFrame)]
    struct F;

    type P = Position<C, F, Meter>;

    fn p(x: f64, y: f64, z: f64) -> P {
        P::new(x, y, z)
    }
    fn m(x: f64) -> Meters {
        Quantity::<Meter>::new(x)
    }

    #[test]
    fn eclipse_umbra_when_fully_covered() {
        let sat = p(0.0, 0.0, 0.0);
        let sun = p(0.0, 0.0, 1.5e11);
        let body = p(0.0, 0.0, 1.0e7);
        assert_eq!(
            eclipse_state(&sat, &sun, m(6.96e8), &body, m(6.378e8)),
            EclipseState::Umbra
        );
    }

    #[test]
    fn eclipse_sunlight_when_unobstructed() {
        let sat = p(0.0, 0.0, 0.0);
        let sun = p(0.0, 0.0, 1.5e11);
        let body = p(0.0, 1.0e10, 0.0);
        assert_eq!(
            eclipse_state(&sat, &sun, m(6.96e8), &body, m(6.378e6)),
            EclipseState::Sunlight
        );
    }
}
