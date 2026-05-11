// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Force models for spacecraft dynamics.
//!
//! ## Trait design
//!
//! [`ForceModel<C, F>`] is the central abstraction.  Each call to
//! [`acceleration`][ForceModel::acceleration] or [`partials`][ForceModel::partials]
//! receives a [`DynamicsContext`] reference so force models can query
//! ephemeris, atmosphere, and gravity-field providers **without holding them
//! as fields**.  Force models store only their tunable physical parameters
//! (C_D, C_R, A/m, J2 coefficient, truncation degree, …).
//!
//! ## Partial derivatives
//!
//! [`ForcePartials<F>`] is the linearization block
//! `A(t) = ∂a/∂[r, v]`, the lower half of the variational matrix
//! `F(t) = [[0, I], [A_r, A_v]]`:
//!
//! * `d_acc_d_pos` (`A_r = ∂a/∂r`): units km/s² per km = s⁻² (stored as raw f64 in frame `F`)
//! * `d_acc_d_vel` (`A_v = ∂a/∂v`): units km/s² per km/s = s⁻¹ (stored as raw f64 in frame `F`)
//!
//! [`TwoBody`] and [`J2`] supply analytic partials.  All other models return
//! `Err(DynamicsError::Provider(_))` from the default impl.
//!
//! ## Provided models
//!
//! | Model | Type parameters | Description |
//! |-------|-----------------|-------------|
//! | [`TwoBody`] | `<Geocentric, GCRS>` | Central Newtonian gravity |
//! | [`J2`] | `<Geocentric, GCRS>` | Zonal oblateness perturbation |
//! | [`DragForce<D>`] | `<Geocentric, GCRS>` | Cannonball atmospheric drag |
//! | [`ThirdBodySunMoon`] | `<Geocentric, GCRS>` | Sun + Moon point-mass perturbation |
//! | [`CannonballSrp`] | `<Geocentric, GCRS>` | Cannonball solar radiation pressure |
//! | [`CompositeForce<C, F>`] | generic | Linear sum of any force models |
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §1, §8, §8.8.1.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2, §3.4.

pub mod composite;
pub mod drag;
pub mod j2;
pub mod srp;
pub mod third_body;
pub mod traits;
pub mod two_body;

// Re-export the full public surface so callers using `use forces::*` or
// `use forces::TwoBody` continue to work unchanged.
pub use composite::CompositeForce;
pub use drag::{DragForce, ExponentialDrag};
pub use j2::J2;
pub use srp::CannonballSrp;
pub use third_body::ThirdBodySunMoon;
pub use traits::{
    ForceModel, ForcePartials, GM_EARTH, GM_MOON, GM_SUN, MU_MOON, MU_SUN, OMEGA_EARTH_RAD_S, P0,
    R_EARTH,
};
pub use two_body::TwoBody;

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
    fn composite_sums_components() {
        let ctx = DynamicsContext::empty();
        let f = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(J2::earth()));
        let a = f.acceleration(&leo(), &ctx).unwrap();
        assert!(a.x().value() < 0.0);
        assert!(a.z().value().abs() < 1e-12);
    }
}
