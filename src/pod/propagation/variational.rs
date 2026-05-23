// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Variational equations for **estimable parameters** beyond the 6×6 state STM.
//!
//! Upstream `siderust::astro::dynamics::variational` integrates the state
//! transition matrix `Φ(t,t₀) = ∂y(t)/∂y(t₀)` with analytic partials. POD
//! also needs the *parameter* partials `S(t) = ∂y(t)/∂p` for each estimable
//! force-model parameter `p` (drag scale Cd_scale, SRP scale Crp_scale,
//! empirical-acceleration coefficients, …).
//!
//! Some parameters do not have closed-form partials of `∂a/∂p` exposed by
//! the upstream force-model trait. The robust, model-agnostic fallback used
//! here is **central finite differences**:
//!
//! ```text
//! ∂y(t)/∂p ≈ [y(t; p + ε) − y(t; p − ε)] / (2 ε)
//! ```
//!
//! The caller supplies a closure that maps a perturbed parameter value to a
//! freshly-built [`ForceModel`]. This keeps the harness completely
//! force-model agnostic — drag scale, SRP scale, and empirical-acceleration
//! amplitudes all reduce to "vary one scalar, rebuild the force, re-propagate".
//!
//! # Validation gate
//!
//! Phase 4 requires that the central-difference column for the drag-scale
//! parameter agree with a *Richardson-refined* estimate to better than
//! `1e-6` relative. [`param_partials_central_diff`] returns the column at
//! step `ε`; the test in `tests/` calls it twice (with `ε` and `ε/2`) and
//! checks that the two estimates agree.

use crate::astro::dynamics::{DynamicsContext, OrbitState, StateTransitionMatrix};
use crate::coordinates::frames::GCRS;
use qtty::Second;

use crate::pod::force::SiderustAccelerationModel;
use crate::pod::propagation::Integrator;
use affn::matrix6::FrameMatrix6;

use super::pod_error::PodDynamicsError;

/// Six-component parameter sensitivity column `∂y/∂p`.
///
/// Layout: `[∂rx/∂p, ∂ry/∂p, ∂rz/∂p, ∂vx/∂p, ∂vy/∂p, ∂vz/∂p]` in km / unit
/// of `p` and km/s / unit of `p`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParamColumn(pub [f64; 6]);

impl ParamColumn {
    /// Element-wise `|self - other|.max()` divided by the max magnitude in
    /// `other`. Returns `0.0` when `other` is identically zero.
    pub fn max_relative_difference(&self, other: &ParamColumn) -> f64 {
        let scale = other.0.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let scale = if scale > 0.0 { scale } else { 1.0 };
        self.0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| (a - b).abs() / scale)
            .fold(0.0_f64, f64::max)
    }
}

/// Diagnostic report bundling a coarse and a refined finite-difference
/// estimate of `∂y/∂p`, plus the relative agreement between them. Used by
/// the validation tests as a Richardson convergence check.
#[derive(Debug, Clone, Copy)]
pub struct ParamStmReport {
    /// Column estimated at step `ε`.
    pub coarse: ParamColumn,
    /// Column estimated at step `ε / 2`.
    pub refined: ParamColumn,
    /// Maximum element-wise relative difference between `coarse` and `refined`.
    pub max_rel_diff: f64,
}

/// Compute the parameter-sensitivity column `∂y(t)/∂p` by central
/// finite differences.
///
/// `build_force(p)` must rebuild the force model from a scalar parameter
/// `p`; the function evaluates it at `p ± ε` and integrates each branch
/// from `state` over `dt` with the supplied [`Integrator`].
///
/// # Arguments
///
/// * `integrator` — fixed-step or adaptive integrator (any [`Integrator`]).
/// * `state` — initial [`OrbitState`] at `t₀`.
/// * `dt` — propagation interval.
/// * `ctx` — dynamics context (atmosphere/ephemeris/gravity providers).
/// * `p_nominal` — nominal value of the parameter.
/// * `epsilon` — perturbation step. Must be strictly positive and finite.
/// * `build_force` — closure `(p) -> impl ForceModel` that builds the
///   force-model composite at parameter value `p`.
///
/// # Errors
///
/// * [`PodDynamicsError::InvalidParamStep`] if `epsilon ≤ 0` or non-finite.
/// * Any [`PodDynamicsError::Dynamics`] raised by the integrator.
///
/// # Example
///
/// ```
/// use siderust::pod::propagation::{param_partials_central_diff, Rk4Integrator};
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::{OrbitState, DynamicsContext};
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::Second;
/// // Trivial parameter: scale GM (no estimable physical meaning, just a smoke test).
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let col = param_partials_central_diff(
///     &Rk4Integrator { step: Second::new(10.0) },
///     s0, Second::new(60.0), &DynamicsContext::empty(),
///     1.0, 1e-6,
///     |_scale| TwoBody::new(siderust::astro::dynamics::GM_EARTH),
/// ).unwrap();
/// assert!(col.0.iter().all(|v| v.is_finite()));
/// ```
pub fn param_partials_central_diff<I, FM, B>(
    integrator: &I,
    state: OrbitState,
    dt: Second,
    ctx: &DynamicsContext,
    p_nominal: f64,
    epsilon: f64,
    mut build_force: B,
) -> Result<ParamColumn, PodDynamicsError>
where
    I: Integrator,
    FM: SiderustAccelerationModel,
    B: FnMut(f64) -> FM,
{
    if !epsilon.is_finite() || epsilon <= 0.0 {
        return Err(PodDynamicsError::InvalidParamStep(
            "epsilon must be strictly positive and finite",
        ));
    }
    let force_plus = build_force(p_nominal + epsilon);
    let force_minus = build_force(p_nominal - epsilon);
    let y_plus = integrator.propagate(&force_plus, state, dt, ctx)?;
    let y_minus = integrator.propagate(&force_minus, state, dt, ctx)?;
    let inv_2eps = 1.0 / (2.0 * epsilon);
    let col = [
        (y_plus.position.x().value() - y_minus.position.x().value()) * inv_2eps,
        (y_plus.position.y().value() - y_minus.position.y().value()) * inv_2eps,
        (y_plus.position.z().value() - y_minus.position.z().value()) * inv_2eps,
        (y_plus.velocity.x().value() - y_minus.velocity.x().value()) * inv_2eps,
        (y_plus.velocity.y().value() - y_minus.velocity.y().value()) * inv_2eps,
        (y_plus.velocity.z().value() - y_minus.velocity.z().value()) * inv_2eps,
    ];
    Ok(ParamColumn(col))
}

// ─────────────────────────────────────────────────────────────────────────────
// PropagatedArc
// ─────────────────────────────────────────────────────────────────────────────

/// A propagated orbit arc with a state-transition matrix accumulated from
/// the initial epoch.
///
/// Each element `(state, Phi)` of [`PropagatedArc::steps`] represents the
/// orbit state and the *cumulative* STM `Φ(tₖ, t₀)` at the k-th step time.
/// Step `0` (index `0`) corresponds to `t₀ + step` with `Φ(t₁, t₀)`.
///
/// # Example
///
/// ```
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::{DynamicsContext, OrbitState, Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use qtty::Second;
/// use siderust::time::JulianDate;
/// use siderust::pod::propagation::{PropagatedArc, VariationalPropagator};
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.545, 0.0),
/// );
/// let prop = VariationalPropagator { step: Second::new(30.0) };
/// let arc = prop.propagate(&TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0, 3, &DynamicsContext::empty()).unwrap();
/// assert_eq!(arc.steps.len(), 3);
/// // Cumulative STM at first step — diagonal must remain close to 1.
/// let phi1 = arc.steps[0].1.as_array();
/// assert!((phi1[0][0] - 1.0).abs() < 0.01);
/// ```
#[derive(Debug, Clone)]
pub struct PropagatedArc {
    /// Ordered list of `(state, Φ(tₖ, t₀))` for each integration step.
    pub steps: Vec<(OrbitState, StateTransitionMatrix<GCRS>)>,
}

impl PropagatedArc {
    /// Number of propagated steps.
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    /// `true` iff no steps were accumulated.
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Final orbit state at the end of the arc (last step), or `None` when
    /// the arc is empty.
    pub fn final_state(&self) -> Option<&OrbitState> {
        self.steps.last().map(|(s, _)| s)
    }

    /// Final cumulative STM `Φ(t_end, t₀)`, or `None` when the arc is empty.
    pub fn final_stm(&self) -> Option<&StateTransitionMatrix<GCRS>> {
        self.steps.last().map(|(_, phi)| phi)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VariationalPropagator
// ─────────────────────────────────────────────────────────────────────────────

/// Propagates an orbit arc step-by-step, accumulating the state-transition
/// matrix `Φ(tₖ, t₀)` at each step.
///
/// Uses [`siderust::astro::dynamics::propagate_stm`] for each sub-step.
/// The cumulative STM is updated by left-multiplying the step STM:
/// `Φ(tₖ₊₁, t₀) = Φ(tₖ₊₁, tₖ) · Φ(tₖ, t₀)`.
///
/// # Example
///
/// ```
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::{DynamicsContext, OrbitState, Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use qtty::Second;
/// use siderust::time::JulianDate;
/// use siderust::pod::propagation::VariationalPropagator;
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.545, 0.0),
/// );
/// let prop = VariationalPropagator { step: Second::new(60.0) };
/// let arc = prop.propagate(&TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0, 5, &DynamicsContext::empty()).unwrap();
/// assert_eq!(arc.len(), 5);
/// assert!(arc.final_state().is_some());
/// ```
pub struct VariationalPropagator {
    /// Fixed time step per integration interval.
    pub step: Second,
}

impl VariationalPropagator {
    /// Propagate `initial` state forward by `n_steps` × `self.step` seconds,
    /// accumulating the cumulative STM at every step.
    ///
    /// The force model must implement the upstream
    /// `siderust::astro::dynamics::forces::ForceModel` trait and must supply
    /// analytic partial derivatives (returned by `ForceModel::partials`) for
    /// the STM to be meaningful. Models without partials produce an identity
    /// STM at every step.
    ///
    /// # Errors
    ///
    /// Returns [`PodDynamicsError::Dynamics`] if `propagate_stm` fails at any step.
    pub fn propagate<FM>(
        &self,
        force: &FM,
        initial: OrbitState,
        n_steps: usize,
        ctx: &DynamicsContext,
    ) -> Result<PropagatedArc, PodDynamicsError>
    where
        FM: SiderustAccelerationModel,
    {
        let mut state = initial;
        let mut phi_total: StateTransitionMatrix<GCRS> = FrameMatrix6::identity();
        let mut steps = Vec::with_capacity(n_steps);

        for _ in 0..n_steps {
            let (new_state, phi_step) = principia::propagate_stm(force, state, self.step, ctx)
                .map_err(|e: principia::PrincipiaError| PodDynamicsError::Dynamics(e.into()))?;
            phi_total = phi_step * phi_total;
            steps.push((new_state, phi_total));
            state = new_state;
        }

        Ok(PropagatedArc { steps })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::pod::propagation::Rk4Integrator;
    use crate::time::JulianDate;

    fn s0() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7_000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
        )
    }

    #[test]
    fn rejects_bad_epsilon() {
        let r = param_partials_central_diff(
            &Rk4Integrator {
                step: Second::new(10.0),
            },
            s0(),
            Second::new(60.0),
            &DynamicsContext::empty(),
            1.0,
            0.0,
            |_| TwoBody::new(crate::astro::dynamics::GM_EARTH),
        );
        assert!(matches!(r, Err(PodDynamicsError::InvalidParamStep(_))));
    }

    #[test]
    fn invariant_force_yields_zero_column() {
        let col = param_partials_central_diff(
            &Rk4Integrator {
                step: Second::new(10.0),
            },
            s0(),
            Second::new(60.0),
            &DynamicsContext::empty(),
            1.0,
            1e-6,
            |_| TwoBody::new(crate::astro::dynamics::GM_EARTH),
        )
        .unwrap();
        assert!(col.0.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn variational_propagator_accumulates_steps() {
        let s = s0();
        let prop = VariationalPropagator {
            step: Second::new(30.0),
        };
        let arc = prop
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s,
                4,
                &DynamicsContext::empty(),
            )
            .unwrap();
        assert_eq!(arc.len(), 4);
        assert!(!arc.is_empty());
        assert!(arc.final_state().is_some());
        assert!(arc.final_stm().is_some());
    }

    #[test]
    fn variational_propagator_stm_close_to_identity_for_short_arc() {
        let s = s0();
        let prop = VariationalPropagator {
            step: Second::new(1.0),
        };
        let arc = prop
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s,
                1,
                &DynamicsContext::empty(),
            )
            .unwrap();
        let phi = arc.steps[0].1.as_array();
        // For dt = 1 s the diagonal of the STM should be very close to 1.
        for (i, row) in phi.iter().enumerate() {
            assert!(
                (row[i] - 1.0).abs() < 1e-2,
                "diagonal entry [{i}][{i}] = {} should be ≈ 1",
                row[i]
            );
        }
    }

    #[test]
    fn prop_arc_empty_when_zero_steps() {
        let s = s0();
        let prop = VariationalPropagator {
            step: Second::new(60.0),
        };
        let arc = prop
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s,
                0,
                &DynamicsContext::empty(),
            )
            .unwrap();
        assert!(arc.is_empty());
        assert!(arc.final_state().is_none());
        assert!(arc.final_stm().is_none());
    }
}
