// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spacecraft state helpers layered on top of `principia`'s Cartesian state.
//!
//! `OrbitState<C, F>` is the astronomy-facing alias for
//! [`principia::DynamicsState<TT, C, F>`]. `siderust` fixes the propagated time
//! scale to [`TT`] while leaving center and frame typed.

use crate::coordinates::cartesian;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::force::Newton;
use crate::qtty::unit::Kilometer;
use crate::qtty::{
    AreaToMass, DragCoefficient, Kilograms, KmPerSecond, KmPerSecondSquared, SquareMeters,
    SrpCoefficient,
};
use crate::time::TT;
use principia::DynamicsState;

/// Geocentric inertial position alias used by most spacecraft APIs.
pub type Position<F = GCRS, U = Kilometer> = cartesian::Position<Geocentric, F, U>;
/// Inertial velocity vector, default `km/s` in [`GCRS`].
pub type Velocity<F = GCRS, U = KmPerSecond> = cartesian::Velocity<F, U>;
/// Inertial acceleration vector, default `km/s²` in [`GCRS`].
pub type Acceleration<F = GCRS, U = KmPerSecondSquared> = affn::cartesian::Acceleration<F, U>;
/// Frame-tagged force vector, default Newton in [`GCRS`].
pub type Force<F = GCRS, U = Newton> = affn::cartesian::Force<F, U>;
/// Default unit of [`Velocity`] used by the propagator (`km/s`).
pub type VelocityUnit = KmPerSecond;
/// Default unit of [`Acceleration`] used by the propagator (`km/s²`).
pub type AccelerationUnit = KmPerSecondSquared;

/// Astronomy-facing alias for the propagated Cartesian state.
pub type OrbitState<C = Geocentric, F = GCRS> = DynamicsState<TT, C, F>;

/// Physical properties of a spacecraft, used by perturbation models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacecraftProperties {
    /// Total wet mass (kg).
    pub mass: Kilograms,
    /// Drag cross-section (m²).
    pub drag_area: SquareMeters,
    /// Drag coefficient C_D.
    pub cd: DragCoefficient,
    /// SRP cross-section (m²).
    pub srp_area: SquareMeters,
    /// SRP coefficient C_R.
    pub cr: SrpCoefficient,
    /// Precomputed drag area-to-mass ratio (m²/kg).
    pub area_to_mass_drag: AreaToMass,
    /// Precomputed SRP area-to-mass ratio (m²/kg).
    pub area_to_mass_srp: AreaToMass,
}

impl SpacecraftProperties {
    /// Construct from all physical parameters.
    #[inline]
    pub fn new(
        mass: Kilograms,
        drag_area: SquareMeters,
        cd: DragCoefficient,
        srp_area: SquareMeters,
        cr: SrpCoefficient,
    ) -> Self {
        let area_to_mass_drag = AreaToMass::new(drag_area.value() / mass.value());
        let area_to_mass_srp = AreaToMass::new(srp_area.value() / mass.value());
        Self {
            mass,
            drag_area,
            cd,
            srp_area,
            cr,
            area_to_mass_drag,
            area_to_mass_srp,
        }
    }

    /// Drag area-to-mass ratio: `drag_area / mass` (m²/kg).
    #[inline]
    pub fn drag_area_to_mass(&self) -> AreaToMass {
        self.area_to_mass_drag
    }

    /// SRP area-to-mass ratio: `srp_area / mass` (m²/kg).
    #[inline]
    pub fn srp_area_to_mass(&self) -> AreaToMass {
        self.area_to_mass_srp
    }

    /// Reasonable demo defaults for a small LEO platform.
    pub fn demo_leo() -> Self {
        Self::new(
            Kilograms::new(500.0),
            SquareMeters::new(2.0),
            DragCoefficient::new(2.2),
            SquareMeters::new(2.0),
            SrpCoefficient::new(1.3),
        )
    }
}

/// Combined spacecraft state used by higher-level estimators and adapters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacecraftState {
    /// Inertial Cartesian orbit state.
    pub orbit: OrbitState,
    /// Spacecraft physical properties.
    pub properties: SpacecraftProperties,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames::GCRS;
    use crate::qtty::{AreaToMass, Kilograms, SquareMeters};
    use crate::time::{JulianDate, JD};

    #[test]
    fn orbit_state_alias_carries_tt_epoch() {
        let state = OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        assert!((state.epoch.to::<JD>().raw().value() - 2_451_545.0).abs() < 1e-9);
    }

    #[test]
    fn spacecraft_properties_precompute_area_to_mass() {
        let props = SpacecraftProperties::new(
            Kilograms::new(500.0),
            SquareMeters::new(2.0),
            DragCoefficient::new(2.2),
            SquareMeters::new(3.0),
            SrpCoefficient::new(1.3),
        );
        assert_eq!(props.drag_area_to_mass(), AreaToMass::new(0.004));
        assert_eq!(props.srp_area_to_mass(), AreaToMass::new(0.006));
    }
}
