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
