// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Newtonian two-body central gravity — `TwoBody` force model.
//!
//! ## Scope
//!
//! Provides [`TwoBody`], implementing the Newtonian central gravitational
//! acceleration and its analytic partial derivatives.
//!
//! ## Equations
//!
//! The central-gravity acceleration is
//!
//! ```text
//! a = −μ r / |r|³
//! ```
//!
//! The analytic Jacobian `∂a/∂r` is
//!
//! ```text
//! A_r = −μ/r³ I + 3μ/r⁵ r rᵀ
//! ```
//!
//! `A_v = 0` (no velocity dependence).
//!
//! ## Units
//!
//! Position km, velocity km/s, acceleration km/s², `A_r` s⁻², `A_v` s⁻¹.
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>` (Earth-centred inertial).
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §1, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::dynamics::units::GM_EARTH;
use crate::coordinates::frames::GCRS;

use super::traits::{ForceModel, ForcePartials, GravitationalParameter, DEGENERATE_RADIUS_KM};

/// Newtonian central-gravity acceleration `−μ r / |r|³`.
///
/// This model provides analytic partial derivatives via the
/// [`ForceModel::partials`] method.  All other force models in this module
/// use the default error-returning implementation.
#[derive(Debug, Clone, Copy)]
pub struct TwoBody {
    /// Gravitational parameter `GM` (km³/s²).
    pub gm: GravitationalParameter,
}

impl TwoBody {
    /// Earth two-body field with EGM2008 GM.
    pub fn earth() -> Self {
        Self { gm: GM_EARTH }
    }
}

impl ForceModel for TwoBody {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let r = s.position.distance().value();
        if r < DEGENERATE_RADIUS_KM {
            return Err(DynamicsError::DegenerateGeometry {
                reason: "TwoBody: radius near zero",
            });
        }
        let r2 = r * r;
        let k = -self.gm.value() / (r2 * r);
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            k * s.position.x().value(),
            k * s.position.y().value(),
            k * s.position.z().value(),
        ))
    }

    fn partials(
        &self,
        s: &OrbitState,
        _ctx: &DynamicsContext,
    ) -> Result<ForcePartials<GCRS>, DynamicsError> {
        let r_vec = [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
        ];
        let r = s.position.distance().value();
        if r < DEGENERATE_RADIUS_KM {
            return Err(DynamicsError::DegenerateGeometry {
                reason: "TwoBody partials: radius near zero",
            });
        }
        Ok(ForcePartials::two_body(self.gm, r_vec))
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

    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn two_body_points_radially_inward() {
        let ctx = DynamicsContext::empty();
        let a = TwoBody::earth().acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() < 0.0);
        assert!(a.y().value().abs() < 1e-12);
        assert!(a.z().value().abs() < 1e-12);
    }

    #[test]
    fn two_body_partials_nonzero() {
        let ctx = DynamicsContext::empty();
        let p = TwoBody::earth().partials(&leo(), &ctx).unwrap();
        let arr = p.d_acc_d_pos.as_array();
        assert!(arr[0][0].abs() > 1e-10 || arr[1][1].abs() > 1e-10);
    }

    /// Verify that `TwoBody::partials` (trait method) agrees with finite
    /// differences using a 5-point stencil.  Relative tolerance: 1 × 10⁻⁷.
    #[test]
    fn two_body_partials_fd_vs_analytic() {
        let ctx = DynamicsContext::empty();
        let model = TwoBody::earth();

        // Representative LEO state: inclined orbit.
        let r0 = [4_500.0_f64, 4_500.0, 3_000.0];
        let state = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r0[0], r0[1], r0[2]),
            Velocity::<GCRS>::new(0.0, 6.0, 4.0),
        );

        let analytic = model.partials(&state, &ctx).unwrap();
        let A = analytic.d_acc_d_pos.as_array();

        // 5-point stencil step: 1 km.
        let h = 1.0_f64; // km

        for j in 0..3 {
            let perturb = |delta: f64| -> [f64; 3] {
                let mut pos = r0;
                pos[j] += delta;
                let s = OrbitState::new_at_jd(
                    JulianDate::new(2_451_545.0),
                    Position::<GCRS>::new(pos[0], pos[1], pos[2]),
                    Velocity::<GCRS>::new(0.0, 6.0, 4.0),
                );
                let a = model.acceleration(&s, &ctx).unwrap();
                [a.x().value(), a.y().value(), a.z().value()]
            };

            let am2 = perturb(-2.0 * h);
            let am1 = perturb(-h);
            let ap1 = perturb(h);
            let ap2 = perturb(2.0 * h);

            for i in 0..3 {
                let fd = (-ap2[i] + 8.0 * ap1[i] - 8.0 * am1[i] + am2[i]) / (12.0 * h);
                let an = A[i][j];
                let denom = an.abs().max(fd.abs()).max(1e-30);
                let rel_err = (fd - an).abs() / denom;
                assert!(
                    rel_err < 1e-7,
                    "TwoBody ∂a[{i}]/∂r[{j}]: analytic={an:.6e}, fd={fd:.6e}, rel={rel_err:.2e}"
                );
            }
        }
    }
}
