// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball atmospheric drag acceleration.

use affn::cartesian::Vector;
use principia::{AccelerationModel, PrincipiaError};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::density::geodetic_altitude;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::ext_qtty::dynamics::InverseSecond;
use crate::qtty::unit::Kilometer;
use crate::qtty::{AreaToMass, DragCoefficient, InverseSeconds, Kilometers, Ratios};
use crate::time::TT;

use super::OMEGA_EARTH_RAD_S;

const WGS84_A_KM: f64 = 6_378.137;
const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// Cannonball atmospheric drag force model.
#[derive(Debug, Clone)]
pub struct DragForce {
    /// Drag coefficient Cd (dimensionless, typically 2.0–2.4 for LEO spacecraft).
    pub cd: DragCoefficient,
    /// Effective area-to-mass ratio (m² kg⁻¹).
    pub area_to_mass: AreaToMass,
    /// Angular velocity vector of Earth's rotation in GCRS (rad s⁻¹, default ≈ 7.292 × 10⁻⁵ k̂).
    pub omega_earth: Vector<GCRS, InverseSecond>,
}

impl DragForce {
    /// Construct a drag model with the standard WGS-84 Earth rotation rate.
    pub fn new(cd: DragCoefficient, area_to_mass: AreaToMass) -> Self {
        Self::with_omega(cd, area_to_mass, OMEGA_EARTH_RAD_S)
    }

    /// Construct a drag model with a custom Earth rotation rate (useful for tests).
    pub fn with_omega(
        cd: DragCoefficient,
        area_to_mass: AreaToMass,
        omega_z: InverseSeconds,
    ) -> Self {
        Self {
            cd,
            area_to_mass,
            omega_earth: Vector::<GCRS, InverseSecond>::new(0.0, 0.0, omega_z.value()),
        }
    }
}

impl AccelerationModel<DynamicsContext, TT, Geocentric, GCRS> for DragForce {
    fn name(&self) -> &'static str {
        "drag"
    }

    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, PrincipiaError> {
        let provider = ctx
            .require_atmosphere()
            .map_err(DynamicsError::into_principia)?;
        let pos_typed: &affn::cartesian::Position<Geocentric, GCRS, Kilometer> = &s.position;
        let h = geodetic_altitude(pos_typed, Kilometers::new(WGS84_A_KM), Ratios::new(WGS84_F));
        if h.value() < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface {
                altitude_km: h.value(),
            }
            .into_principia());
        }
        let rho = provider
            .density(h)
            .map_err(DynamicsError::into_principia)?
            .value();
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let oz = self.omega_earth.z().value();
        let omega_cross_r_x = -oz * ry;
        let omega_cross_r_y = oz * rx;
        let vx = s.velocity.x().value() - omega_cross_r_x;
        let vy = s.velocity.y().value() - omega_cross_r_y;
        let vz = s.velocity.z().value();
        let v_mag_m_s = (vx * vx + vy * vy + vz * vz).sqrt() * 1_000.0;
        let pre = -0.5 * self.cd.value() * self.area_to_mass.value() * rho * v_mag_m_s;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            pre * vx,
            pre * vy,
            pre * vz,
        ))
    }
}

/// Type alias: [`DragForce`] with the default exponential atmosphere model.
///
/// Kept for API compatibility; prefer [`DragForce`] directly.
pub type ExponentialDrag = DragForce;

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::density::{
        ConstantDensity, DensityProvider, ExponentialAtmosphere,
    };
    use crate::astro::dynamics::forces::R_EARTH;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::qtty::{AreaToMass, DragCoefficient, KilogramsPerCubicMeter};
    use crate::time::JulianDate;

    fn ctx_with_exponential() -> DynamicsContext {
        let atm: Arc<dyn DensityProvider + Send + Sync> =
            Arc::new(ExponentialAtmosphere::LEO_500KM);
        DynamicsContextBuilder::new().with_atmosphere(atm).build()
    }

    fn ctx_with_constant(rho: f64) -> DynamicsContext {
        let atm: Arc<dyn DensityProvider + Send + Sync> = Arc::new(ConstantDensity {
            rho: KilogramsPerCubicMeter::new(rho),
        });
        DynamicsContextBuilder::new().with_atmosphere(atm).build()
    }

    fn state(radius_km: f64, vy_km_s: f64) -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(radius_km, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, vy_km_s, 0.0),
        )
    }

    #[test]
    fn drag_error_when_no_provider() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let err = drag.acceleration(
            &state(R_EARTH.value() + 400.0, 7.7),
            &DynamicsContext::empty(),
        );
        assert!(matches!(
            err,
            Err(PrincipiaError::ContextDataUnavailable { what: "atmosphere" })
        ));
    }

    #[test]
    fn drag_direction_opposes_velocity() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let acc = drag
            .acceleration(
                &state(R_EARTH.value() + 400.0, 7.7),
                &ctx_with_exponential(),
            )
            .unwrap();
        assert!(acc.y().value() < 0.0);
    }

    #[test]
    fn drag_below_surface_returns_error() {
        let drag = DragForce::new(DragCoefficient::new(2.2), AreaToMass::new(0.01));
        let err = drag.acceleration(
            &state(R_EARTH.value() - 100.0, 7.0),
            &ctx_with_constant(1.0e-9),
        );
        assert!(matches!(
            err,
            Err(PrincipiaError::DegenerateGeometry {
                reason: "altitude below surface"
            })
        ));
    }
}
