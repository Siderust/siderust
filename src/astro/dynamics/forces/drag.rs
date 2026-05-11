// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball atmospheric drag — `DragForce<D>` force model.
//!
//! ## Scope
//!
//! Provides [`DragForce<D>`], a cannonball aerodynamic drag model using an
//! embedded [`DensityProvider`].
//!
//! ## Equations
//!
//! ```text
//! a_drag = − ½ · Cd · (A/m) · ρ(h) · |v_rel| · v_rel
//! ```
//!
//! where `v_rel = v − ω_⊕ × r` is the velocity relative to the co-rotating
//! atmosphere.
//!
//! ## Units
//!
//! Position km, velocity km/s, acceleration km/s², density kg/m³.
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>`.  Earth rotation is modelled via a fixed sidereal
//! rate [`OMEGA_EARTH_RAD_S`].
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.4.

use crate::astro::dynamics::atmosphere::{DensityProvider, ExponentialAtmosphere};
use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;
use crate::qtty::{AreaToMass, DragCoefficient, Kilometers};

use super::traits::{ForceModel, OMEGA_EARTH_RAD_S, R_EARTH};

/// Cannonball atmospheric drag acceleration.
///
/// ```text
/// a_drag = − ½ · Cd · (A/m) · ρ(h) · |v_rel| · v_rel
/// ```
///
/// where:
///
/// * `Cd` is the drag coefficient (typed [`DragCoefficient`]; typical LEO value ≈ 2.2);
/// * `A/m` is the area-to-mass ratio (typed [`AreaToMass`], m²/kg);
/// * `ρ(h)` is the atmospheric mass density from the embedded [`DensityProvider`];
/// * `h = |r| − R_⊕` is the geocentric altitude (km);
/// * `v_rel = v − ω_⊕ × r` is the velocity in the co-rotating atmosphere frame.
///
/// The embedded atmosphere `D` holds the density-model parameters (e.g.
/// exponential-atmosphere scale heights).
#[derive(Debug, Clone)]
pub struct DragForce<D: DensityProvider> {
    /// Drag coefficient C_D (dimensionless).
    pub cd: DragCoefficient,
    /// Effective area-to-mass ratio (m²/kg).
    pub area_to_mass: AreaToMass,
    /// Atmosphere density provider.
    pub atmosphere: D,
}

impl DragForce<ExponentialAtmosphere> {
    /// Build a drag model using the [`ExponentialAtmosphere::LEO_500KM`] profile.
    pub fn leo_500km(cd: f64, area_to_mass_m2_kg: f64) -> Self {
        Self {
            cd: DragCoefficient::new(cd),
            area_to_mass: AreaToMass::new(area_to_mass_m2_kg),
            atmosphere: ExponentialAtmosphere::LEO_500KM,
        }
    }
}

impl<D: DensityProvider + Send + Sync> ForceModel for DragForce<D> {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        let h = r - R_EARTH.value();
        if h < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface { altitude_km: h });
        }
        let rho = self.atmosphere.density(Kilometers::new(h)).value();

        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let vx = s.velocity.x().value();
        let vy = s.velocity.y().value();
        let vz = s.velocity.z().value();

        let omega_cross_r = [-OMEGA_EARTH_RAD_S * ry, OMEGA_EARTH_RAD_S * rx, 0.0_f64];
        let v_rel = [
            vx - omega_cross_r[0],
            vy - omega_cross_r[1],
            vz - omega_cross_r[2],
        ];

        let v_mag_m_s = (v_rel[0].powi(2) + v_rel[1].powi(2) + v_rel[2].powi(2)).sqrt() * 1_000.0;
        let pre = -0.5 * self.cd.value() * self.area_to_mass.value() * rho * v_mag_m_s;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            pre * v_rel[0],
            pre * v_rel[1],
            pre * v_rel[2],
        ))
    }
}

/// Type alias: drag model with the built-in exponential atmosphere.
pub type ExponentialDrag = DragForce<ExponentialAtmosphere>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::atmosphere::ExponentialAtmosphere;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::integrators::rk4_propagate;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::qtty::{AreaToMass, DragCoefficient, KilogramsPerCubicMeter, Kilometers, Second};
    use crate::time::JulianDate;

    use super::super::composite::CompositeForce;
    use super::super::traits::R_EARTH;
    use super::super::two_body::TwoBody;

    #[test]
    fn drag_density_decreases_with_altitude() {
        let d = DragForce::leo_500km(2.2, 0.02);
        assert!(
            d.atmosphere.density(Kilometers::new(500.0))
                > d.atmosphere.density(Kilometers::new(600.0))
        );
        assert!(
            d.atmosphere.density(Kilometers::new(400.0))
                > d.atmosphere.density(Kilometers::new(500.0))
        );
    }

    #[test]
    fn drag_decays_orbit_altitude() {
        let mu: f64 = 398_600.441_8;
        let r0 = R_EARTH.value() + 350.0;
        let v0 = (mu / r0).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v0, 0.0),
        );
        let force = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(DragForce {
                cd: DragCoefficient::new(2.2),
                area_to_mass: AreaToMass::new(5.0),
                atmosphere: ExponentialAtmosphere {
                    rho0: KilogramsPerCubicMeter::new(1.0e-11),
                    h0: Kilometers::new(350.0),
                    scale_height: Kilometers::new(50.0),
                },
            }));
        let ctx = DynamicsContext::empty();
        let s_end = rk4_propagate(&force, s0, Second::new(30.0), 360, &ctx).unwrap();
        let r_end = s_end.position.distance().value();
        assert!(
            r_end < r0,
            "expected drag-driven decay; r0={r0:.3}, r_end={r_end:.3}"
        );
    }
}
