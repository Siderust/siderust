// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! `CompositeForce` — linear sum of N component force models.
//!
//! ## Scope
//!
//! Provides [`CompositeForce<C, F>`], which accumulates the accelerations and
//! partial derivatives of an ordered list of [`ForceModel`] implementations.
//!
//! ## Equations
//!
//! For components `{a₁, a₂, …, aₙ}`:
//!
//! ```text
//! a_total   = Σ aᵢ
//! A_r_total = Σ ∂aᵢ/∂r
//! A_v_total = Σ ∂aᵢ/∂v
//! ```
//!
//! ## Units
//!
//! km/s² (acceleration), s⁻² (`A_r`), s⁻¹ (`A_v`).
//!
//! ## Frame/center assumptions
//!
//! Generic over `<C, F>`; defaults to `<Geocentric, GCRS>`.
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};

use super::traits::{ForceModel, ForcePartials};

/// Sum of N component force models.
///
/// `acceleration` sums child accelerations, short-circuiting on the first
/// error.  `partials` sums child [`ForcePartials`], also short-circuiting on
/// the first error.
///
/// # Type parameters
///
/// * `C` — reference center (default [`Geocentric`]).
/// * `F` — reference frame (default [`GCRS`]).
pub struct CompositeForce<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    components: Vec<Box<dyn ForceModel<C, F>>>,
}

impl<C, F> Default for CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    fn default() -> Self {
        Self::empty()
    }
}

impl<C, F> CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    /// Empty composite force (returns zero acceleration).
    pub fn empty() -> Self {
        Self {
            components: Vec::new(),
        }
    }

    /// Append a force component.
    pub fn push(mut self, f: Box<dyn ForceModel<C, F>>) -> Self {
        self.components.push(f);
        self
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// True when no components have been added.
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
}

impl<C, F> ForceModel<C, F> for CompositeForce<C, F>
where
    C: ReferenceCenter + 'static,
    F: ReferenceFrame + 'static,
{
    fn acceleration(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<F, AccelerationUnit>, DynamicsError> {
        let mut acc = Acceleration::<F, AccelerationUnit>::new(0.0, 0.0, 0.0);
        for f in &self.components {
            acc = acc + f.acceleration(state, ctx)?;
        }
        Ok(acc)
    }

    fn partials(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<ForcePartials<F>, DynamicsError> {
        let mut total = ForcePartials::zero();
        for f in &self.components {
            total.add_in_place(&f.partials(state, ctx)?);
        }
        Ok(total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::traits::ForceModel;
    use crate::astro::dynamics::forces::two_body::TwoBody;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn sample_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn empty_composite_is_zero() {
        let comp = CompositeForce::empty();
        assert!(comp.is_empty());
        assert_eq!(comp.len(), 0);
        let ctx = DynamicsContext::empty();
        let a = comp.acceleration(&sample_state(), &ctx).unwrap();
        assert_eq!(a.x().value(), 0.0);
        assert_eq!(a.y().value(), 0.0);
        assert_eq!(a.z().value(), 0.0);
    }

    #[test]
    fn default_composite_is_empty() {
        let comp: CompositeForce = CompositeForce::default();
        assert!(comp.is_empty());
    }

    #[test]
    fn push_increases_len() {
        let comp = CompositeForce::empty().push(Box::new(TwoBody::earth()));
        assert_eq!(comp.len(), 1);
        assert!(!comp.is_empty());
    }

    #[test]
    fn composite_acceleration_sums_components() {
        let comp = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(TwoBody::earth()));
        let ctx = DynamicsContext::empty();
        let a = comp.acceleration(&sample_state(), &ctx).unwrap();
        let single = TwoBody::earth().acceleration(&sample_state(), &ctx).unwrap();
        let rel_x = (a.x().value() - 2.0 * single.x().value()).abs() / single.x().value().abs();
        assert!(rel_x < 1e-12, "composite x should be 2× single-component: got {rel_x}");
    }

    #[test]
    fn composite_partials_sums_components() {
        let comp = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(TwoBody::earth()));
        let ctx = DynamicsContext::empty();
        let p = comp.partials(&sample_state(), &ctx).unwrap();
        let single = TwoBody::earth().partials(&sample_state(), &ctx).unwrap();
        let pa = p.d_acc_d_pos.as_array();
        let sa = single.d_acc_d_pos.as_array();
        for i in 0..3 {
            for j in 0..3 {
                let rel = (pa[i][j] - 2.0 * sa[i][j]).abs();
                assert!(rel < 1e-12, "partials[{i}][{j}] mismatch: {rel}");
            }
        }
    }
}
