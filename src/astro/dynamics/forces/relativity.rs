// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! First-order post-Newtonian (Schwarzschild) relativistic correction.
//!
//! ## Scope
//!
//! Provides [`CentralBodyRelativity1Pn`], implementing the Schwarzschild
//! relativistic acceleration for a central body.
//!
//! ## Equation
//!
//! IERS 2010 §10.3, Schwarzschild term (geocentric, neglecting J2 and
//! Lense–Thirring):
//!
//! ```text
//! a_rel = (GM / (c² r³)) [ (4GM/r − v·v) r + 4(r·v) v ]
//! ```
//!
//! where `r` and `v` are the position (km) and velocity (km/s) vectors and
//! `c` is the speed of light (km/s).
//!
//! ## Units
//!
//! Position km, velocity km/s, acceleration km/s².
//!
//! ## References
//!
//! * IERS Conventions (2010), §10.3.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::frames::GCRS;
use crate::qtty::{KmPerSeconds, SPEED_OF_LIGHT_KM_S};

use super::traits::{ForceModel, GravitationalParameter, GM_EARTH};

/// First-order post-Newtonian (Schwarzschild) relativistic correction for a
/// central body.
///
/// This model adds the general-relativistic precession-inducing correction to
/// the two-body acceleration.  The Lense–Thirring (frame-dragging) and J2
/// coupling terms are neglected.
#[derive(Debug, Clone, Copy)]
pub struct CentralBodyRelativity1Pn {
    /// Central-body gravitational parameter `GM` (km³/s²).
    pub gm: GravitationalParameter,
    /// Speed of light (km/s).
    pub c: KmPerSeconds,
}

impl CentralBodyRelativity1Pn {
    /// Earth Schwarzschild correction using EGM2008 GM and the exact SI speed
    /// of light.
    pub fn earth() -> Self {
        Self {
            gm: GM_EARTH,
            c: SPEED_OF_LIGHT_KM_S,
        }
    }
}

impl ForceModel for CentralBodyRelativity1Pn {
    /// Schwarzschild relativistic acceleration (IERS 2010 §10.3).
    ///
    /// ```text
    /// a_rel = (GM / (c² r³)) [ (4GM/r − v·v) r + 4(r·v) v ]
    /// ```
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let rx = s.position.x().value();
        let ry = s.position.y().value();
        let rz = s.position.z().value();

        let vx = s.velocity.x().value();
        let vy = s.velocity.y().value();
        let vz = s.velocity.z().value();

        let r2 = rx * rx + ry * ry + rz * rz;
        let r = r2.sqrt();
        let r3 = r2 * r;

        let v2 = vx * vx + vy * vy + vz * vz;
        let rdotv = rx * vx + ry * vy + rz * vz;

        let gm = self.gm.value();
        let c = self.c.value();
        let c2 = c * c;

        // scalar prefactor: GM / (c² r³)
        let pre = gm / (c2 * r3);

        // bracket: (4GM/r - v·v)
        let coeff_r = 4.0 * gm / r - v2;
        // bracket: 4(r·v)
        let coeff_v = 4.0 * rdotv;

        let ax = pre * (coeff_r * rx + coeff_v * vx);
        let ay = pre * (coeff_r * ry + coeff_v * vy);
        let az = pre * (coeff_r * rz + coeff_v * vz);

        Ok(Acceleration::<GCRS, AccelerationUnit>::new(ax, ay, az))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn leo_circular() -> OrbitState {
        // Circular LEO: r = 7000 km, v ≈ 7.5 km/s along +Y, position along +X.
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn relativity_magnitude_order_of_magnitude() {
        let ctx = DynamicsContext::empty();
        let a = CentralBodyRelativity1Pn::earth()
            .acceleration(&leo_circular(), &ctx)
            .unwrap();
        let mag = (a.x().value().powi(2) + a.y().value().powi(2) + a.z().value().powi(2)).sqrt();
        assert!(
            mag > 1e-12 && mag < 1e-9,
            "expected 1e-12 < |a_rel| < 1e-9 km/s², got {mag:.3e}"
        );
    }

    #[test]
    fn relativity_approximately_radial_for_circular_orbit() {
        let ctx = DynamicsContext::empty();
        let state = leo_circular();
        let a = CentralBodyRelativity1Pn::earth()
            .acceleration(&state, &ctx)
            .unwrap();

        // Radial unit vector
        let rx = state.position.x().value();
        let ry = state.position.y().value();
        let rz = state.position.z().value();
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        let r_hat = [rx / r, ry / r, rz / r];

        let ax = a.x().value();
        let ay = a.y().value();
        let az = a.z().value();
        let a_mag = (ax * ax + ay * ay + az * az).sqrt();

        // cosine of angle between a and r̂
        let cos_theta = (ax * r_hat[0] + ay * r_hat[1] + az * r_hat[2]) / a_mag;
        assert!(
            cos_theta.abs() > 0.9,
            "relativity acceleration should be roughly radial for a circular orbit, cos={cos_theta:.4}"
        );
    }
}
