// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball solar radiation pressure — `CannonballSrp` force model.
//!
//! ## Scope
//!
//! Provides [`CannonballSrp`], a cannonball SRP model that always treats the
//! spacecraft as fully sunlit (no eclipse modelling yet).
//!
//! ## Equations
//!
//! ```text
//! a_srp = Cr · P₀ · (AU / |r_sun_sat|)² · (A/m) · r̂_sun_sat
//! ```
//!
//! The acceleration pushes the spacecraft away from the Sun.
//!
//! ## Units
//!
//! Position km, P₀ N/m², area-to-mass m²/kg → acceleration km/s² after
//! dividing by 1 000.
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>`.
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.4.

use affn::cartesian::Displacement;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;
use crate::qtty::{AreaToMass, Kilometer, SrpCoefficient};

use super::third_body::sun_geocentric;
use super::traits::{ForceModel, AU_IN_KM, P0};

/// Cannonball solar radiation pressure (SRP).
///
/// ```text
/// a_srp = Cr · P0 · (AU / |r_sun_sat|)² · (A/m) · r̂_sun_sat
/// ```
///
/// The resulting acceleration pushes the spacecraft away from the Sun.
/// Eclipse modelling is not yet included — the Sun is always treated as visible.
/// The Sun geocentric position is fetched from the [`DynamicsContext`] ephemeris.
#[derive(Debug, Clone, Copy)]
pub struct CannonballSrp {
    /// Radiation-pressure coefficient (dimensionless).
    pub cr: SrpCoefficient,
    /// Area-to-mass ratio (m²/kg).
    pub area_to_mass: AreaToMass,
}

impl CannonballSrp {
    /// Build a cannonball SRP model.
    pub fn new(cr: f64, area_to_mass_m2_kg: f64) -> Self {
        Self {
            cr: SrpCoefficient::new(cr),
            area_to_mass: AreaToMass::new(area_to_mass_m2_kg),
        }
    }
}

impl ForceModel for CannonballSrp {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let eph = ctx.require_ephemeris()?;
        let sun = sun_geocentric(eph, s.epoch_jd())?;
        let r_sun_sat = Displacement::<GCRS, Kilometer>::new(
            s.position.x().value() - sun.x().value(),
            s.position.y().value() - sun.y().value(),
            s.position.z().value() - sun.z().value(),
        );
        let r = r_sun_sat.magnitude().value();
        if r == 0.0 {
            return Ok(Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0));
        }
        let r2 = r * r;
        // N/m² · m²/kg = m/s²; convert to km/s² by dividing by 1000.
        let mag_km_s2 =
            self.cr.value() * P0.value() * (AU_IN_KM * AU_IN_KM / r2) * self.area_to_mass.value()
                / 1_000.0;
        let inv_r = 1.0 / r;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            mag_km_s2 * r_sun_sat.x().value() * inv_r,
            mag_km_s2 * r_sun_sat.y().value() * inv_r,
            mag_km_s2 * r_sun_sat.z().value() * inv_r,
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::calculus::ephemeris::Vsop87Ephemeris;
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn leo_at(epoch: JulianDate) -> OrbitState {
        OrbitState::new_at_jd(
            epoch,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn srp_order_of_magnitude_at_leo() {
        let ctx = DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build();
        let srp = CannonballSrp::new(1.5, 0.02);
        let s = leo_at(JulianDate::new(2_451_545.0));
        let a = srp.acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(
            (5e-11..5e-10).contains(&mag),
            "SRP magnitude out of expected band: {mag} km/s²"
        );
    }

    #[test]
    fn srp_zero_when_area_is_zero() {
        let ctx = DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build();
        let srp = CannonballSrp::new(1.5, 0.0);
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() == 0.0 && a.y().value() == 0.0 && a.z().value() == 0.0);
    }
}
