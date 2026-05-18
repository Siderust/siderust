// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Variational-equations propagator — analytic STM via RK4.
//!
//! ## Entry points
//!
//! - [`propagate_stm`] — propagate `(orbit state, STM)` over a signed interval
//!   using an automatically-chosen fixed step.
//! - [`propagate_stm_with`] — same, with an explicit [`VariationalConfig`].
//!
//! ## Algorithm
//!
//! Both functions use a classical RK4 scheme applied jointly to the 6-DOF
//! orbit state and the 6×6 STM (42 scalars in total). Each sub-step:
//!
//! 1. Calls `ForceModel::acceleration` for the orbit ODE.
//! 2. Calls `ForceModel::partials` to build the 6×6 Jacobian `A(t)`.
//! 3. Updates Φ via `dΦ/dt = A · Φ`.
//!
//! The orbit state and STM are advanced together through all four RK4 stages,
//! ensuring that intermediate Φ values (used in stages k2, k3, k4) are
//! self-consistent with the corresponding intermediate orbit state.
//!
//! ## Step size
//!
//! The default step is `min(|dt| / 100, 30 s)`. This is conservative and
//! intentionally not adaptive — correctness over performance for this release.
//! Use [`propagate_stm_with`] to override.
//!
//! ## Signed propagation
//!
//! `dt` may be negative (backward propagation). The number of sub-steps is
//! `ceil(|dt| / step)` and each internal step has the same sign as `dt`.

use affn::matrix6::FrameMatrix6;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::traits::ForceModel;
use crate::astro::dynamics::state::{OrbitState, Position, Velocity};
use crate::coordinates::frames::GCRS;
use crate::qtty::Second;

use super::equations::{identity_6x6, mat6_add, mat6_scale, variational_derivative};

// =============================================================================
// Public type alias
// =============================================================================

/// Frame-tagged 6×6 state-transition matrix Φ(t, t₀) = `∂x(t)/∂x(t₀)`.
///
/// This is an alias for [`FrameMatrix6<F>`] with the STM semantics documented
/// in [`super`].
pub type StateTransitionMatrix<F> = FrameMatrix6<F>;

// =============================================================================
// Configuration
// =============================================================================

/// Configuration for the variational-equations propagator.
#[derive(Debug, Clone, Copy)]
pub struct VariationalConfig {
    /// Fixed sub-step size. Must be positive.
    ///
    /// The actual internal step will have the same sign as `dt` so that
    /// forward and backward propagation both work.
    pub step: Second,
}

impl Default for VariationalConfig {
    /// Default step: 30 s (conservative; override for performance-critical use).
    fn default() -> Self {
        Self {
            step: Second::new(30.0),
        }
    }
}

// =============================================================================
// Public API
// =============================================================================

/// Propagate the orbit state and STM over `dt` seconds using a default step.
///
/// The default step size is `min(|dt| / 100, 30 s)`.
///
/// # Errors
///
/// Returns a [`DynamicsError`] if any force-model call fails. In particular,
/// force models that do not implement analytic partials will cause an error;
/// use only models that override [`ForceModel::partials`] (e.g. [`TwoBody`],
/// [`J2`], [`CompositeForce`] of such models).
///
/// [`TwoBody`]: crate::astro::dynamics::forces::TwoBody
/// [`J2`]: crate::astro::dynamics::forces::J2
/// [`CompositeForce`]: crate::astro::dynamics::forces::CompositeForce
pub fn propagate_stm<FM>(
    force: &FM,
    state: OrbitState,
    dt: Second,
    ctx: &DynamicsContext,
) -> Result<(OrbitState, StateTransitionMatrix<GCRS>), DynamicsError>
where
    FM: ForceModel,
{
    let dt_s = dt.value().abs();
    let step_s = (dt_s / 100.0).clamp(1e-12, 30.0);
    let config = VariationalConfig {
        step: Second::new(step_s),
    };
    propagate_stm_with(force, state, dt, ctx, &config)
}

/// Propagate the orbit state and STM with an explicit step configuration.
///
/// See [`propagate_stm`] for details.
pub fn propagate_stm_with<FM>(
    force: &FM,
    state: OrbitState,
    dt: Second,
    ctx: &DynamicsContext,
    config: &VariationalConfig,
) -> Result<(OrbitState, StateTransitionMatrix<GCRS>), DynamicsError>
where
    FM: ForceModel,
{
    let dt_s = dt.value();

    // Zero-interval: return identity STM immediately.
    if dt_s == 0.0 {
        return Ok((state, FrameMatrix6::identity()));
    }

    let step_abs = config.step.value().abs().max(1e-15);
    let n_steps = (dt_s.abs() / step_abs).ceil() as usize;
    let n_steps = n_steps.max(1);
    let h = dt_s / n_steps as f64; // carries the sign of dt

    let mut current_state = state;
    let mut phi = identity_6x6();

    for _ in 0..n_steps {
        let (new_state, new_phi) =
            variational_rk4_step(force, &current_state, &phi, Second::new(h), ctx)?;
        current_state = new_state;
        phi = new_phi;
    }

    Ok((current_state, FrameMatrix6::from_array(phi)))
}

// =============================================================================
// Internal: single RK4 step for the joint (orbit, Φ) system
// =============================================================================

/// One RK4 step for the combined `(orbit_state, Φ)` variational system.
///
/// Both the orbit state and the STM are advanced together so that the
/// intermediate Φ values used in stages k2, k3, k4 remain consistent with
/// the corresponding intermediate orbit states.
fn variational_rk4_step<FM>(
    force: &FM,
    state: &OrbitState,
    phi: &[[f64; 6]; 6],
    h: Second,
    ctx: &DynamicsContext,
) -> Result<(OrbitState, [[f64; 6]; 6]), DynamicsError>
where
    FM: ForceModel,
{
    let h_s = h.value();
    let half = h_s * 0.5;

    // --- Stage 1 ---
    let (dy1, dphi1) = variational_derivative(force, state, phi, ctx)?;

    // --- Stage 2 ---
    let state2 = advance_orbit(state, &dy1, half);
    let phi2 = mat6_add(phi, &mat6_scale(&dphi1, half));
    let (dy2, dphi2) = variational_derivative(force, &state2, &phi2, ctx)?;

    // --- Stage 3 ---
    let state3 = advance_orbit(state, &dy2, half);
    let phi3 = mat6_add(phi, &mat6_scale(&dphi2, half));
    let (dy3, dphi3) = variational_derivative(force, &state3, &phi3, ctx)?;

    // --- Stage 4 ---
    let state4 = advance_orbit(state, &dy3, h_s);
    let phi4 = mat6_add(phi, &mat6_scale(&dphi3, h_s));
    let (dy4, dphi4) = variational_derivative(force, &state4, &phi4, ctx)?;

    // --- Combine: phi_new = phi + h * (k1 + 2k2 + 2k3 + k4) / 6 ---
    let dy_comb = rk4_combine_6(&dy1, &dy2, &dy3, &dy4);
    let dphi_comb = rk4_combine_mat6x6(&dphi1, &dphi2, &dphi3, &dphi4);

    let new_state = advance_orbit_with_epoch(state, &dy_comb, h_s);
    let new_phi = mat6_add(phi, &mat6_scale(&dphi_comb, h_s));

    Ok((new_state, new_phi))
}

// =============================================================================
// Helpers
// =============================================================================

/// Apply Euler step to orbit state (epoch unchanged — used for k2/k3/k4 stages).
#[inline]
fn advance_orbit(state: &OrbitState, dy: &[f64; 6], h: f64) -> OrbitState {
    OrbitState::new(
        state.epoch,
        Position::<GCRS>::new(
            state.position.x().value() + h * dy[0],
            state.position.y().value() + h * dy[1],
            state.position.z().value() + h * dy[2],
        ),
        Velocity::<GCRS>::new(
            state.velocity.x().value() + h * dy[3],
            state.velocity.y().value() + h * dy[4],
            state.velocity.z().value() + h * dy[5],
        ),
    )
}

/// Apply Euler step to orbit state, also advancing the epoch by `h` seconds.
#[inline]
fn advance_orbit_with_epoch(state: &OrbitState, dy: &[f64; 6], h: f64) -> OrbitState {
    let new_epoch = state.epoch + Second::new(h);
    OrbitState::new(
        new_epoch,
        Position::<GCRS>::new(
            state.position.x().value() + h * dy[0],
            state.position.y().value() + h * dy[1],
            state.position.z().value() + h * dy[2],
        ),
        Velocity::<GCRS>::new(
            state.velocity.x().value() + h * dy[3],
            state.velocity.y().value() + h * dy[4],
            state.velocity.z().value() + h * dy[5],
        ),
    )
}

/// RK4 weighted combination for a 6-vector: `(k1 + 2k2 + 2k3 + k4) / 6`.
#[inline]
fn rk4_combine_6(k1: &[f64; 6], k2: &[f64; 6], k3: &[f64; 6], k4: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for i in 0..6 {
        out[i] = (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    out
}

/// RK4 weighted combination for a 6×6 matrix: `(k1 + 2k2 + 2k3 + k4) / 6`.
#[inline]
fn rk4_combine_mat6x6(
    k1: &[[f64; 6]; 6],
    k2: &[[f64; 6]; 6],
    k3: &[[f64; 6]; 6],
    k4: &[[f64; 6]; 6],
) -> [[f64; 6]; 6] {
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = (k1[i][j] + 2.0 * k2[i][j] + 2.0 * k3[i][j] + k4[i][j]) / 6.0;
        }
    }
    out
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::{CompositeForce, TwoBody, J2};
    use crate::astro::dynamics::integrators::rk4::rk4_propagate;
    use crate::astro::dynamics::state::{OrbitState, Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    /// Typical LEO state used across tests.
    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(6_800.0, 1_200.0, 2_400.0),
            Velocity::<GCRS>::new(-1.5, 7.2, 0.8),
        )
    }

    /// Extract the 6-vector `[rx, ry, rz, vx, vy, vz]` from an orbit state.
    fn state_vec(s: &OrbitState) -> [f64; 6] {
        [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
            s.velocity.x().value(),
            s.velocity.y().value(),
            s.velocity.z().value(),
        ]
    }

    /// Perturb an orbit state by a raw 6-vector `δy0`.
    fn perturb(s: &OrbitState, dy: &[f64; 6]) -> OrbitState {
        OrbitState::new(
            s.epoch,
            Position::<GCRS>::new(
                s.position.x().value() + dy[0],
                s.position.y().value() + dy[1],
                s.position.z().value() + dy[2],
            ),
            Velocity::<GCRS>::new(
                s.velocity.x().value() + dy[3],
                s.velocity.y().value() + dy[4],
                s.velocity.z().value() + dy[5],
            ),
        )
    }

    /// Apply Φ to a 6-vector: `Φ · δy0`.
    fn phi_apply(phi: &StateTransitionMatrix<GCRS>, dy: &[f64; 6]) -> [f64; 6] {
        let m = phi.as_array();
        let mut out = [0.0_f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                out[i] += m[i][j] * dy[j];
            }
        }
        out
    }

    /// Compute the relative error between two 6-vectors, normalized by the
    /// infinity norm of `actual` (floored at 1e-30 to avoid zero division).
    ///
    /// This is more meaningful than per-component relative error when some
    /// components are near zero: a picometer error on a 50 μm component is
    /// not a useful signal.
    fn max_rel_err(actual: &[f64; 6], predicted: &[f64; 6]) -> f64 {
        let scale = actual
            .iter()
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max)
            .max(1e-30);
        let mut max = 0.0_f64;
        for i in 0..6 {
            let rel = (actual[i] - predicted[i]).abs() / scale;
            max = max.max(rel);
        }
        max
    }

    /// Determinant of a 6×6 matrix via Gaussian elimination with partial pivoting.
    fn det6(m: [[f64; 6]; 6]) -> f64 {
        let mut a = m;
        let mut det = 1.0_f64;
        for col in 0..6 {
            // Partial pivot
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in (col + 1)..6 {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            if max_row != col {
                a.swap(max_row, col);
                det = -det;
            }
            if a[col][col].abs() < 1e-15 {
                return 0.0;
            }
            det *= a[col][col];
            for row in (col + 1)..6 {
                let factor = a[row][col] / a[col][col];
                for k in col..6 {
                    let sub = factor * a[col][k];
                    a[row][k] -= sub;
                }
            }
        }
        det
    }

    // -------------------------------------------------------------------------

    /// Φ at zero time must be the 6×6 identity within 1e-12 per element.
    #[test]
    fn identity_at_zero_time() {
        let ctx = DynamicsContext::empty();
        let force = TwoBody::earth();
        let s0 = leo();

        let (_, phi) = propagate_stm(&force, s0, Second::new(0.0), &ctx).unwrap();
        let eye = *phi.as_array();
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (eye[i][j] - expected).abs() < 1e-12,
                    "Φ[{i},{j}] = {:.4e}, expected {expected}",
                    eye[i][j]
                );
            }
        }
    }

    /// Two-body STM: Φ · δy0 predicts the propagated displacement to ~1e-6 relative.
    ///
    /// Protocol: perturb initial state by δy0, propagate both nominal and perturbed
    /// states forward 600 s with `rk4_propagate`, then compare
    /// `y(600) - y_nom(600)` to `Φ · δy0` (Φ comes from `propagate_stm`).
    ///
    /// Using the same integrator for both trajectories eliminates floating-point
    /// path differences from the FD reference.  δy0 = 1e-3 km so the second-order
    /// linearization error scales as δy0/r ~ 1e-7 < 1e-6.
    #[test]
    fn two_body_stm_agrees_with_fd() {
        let ctx = DynamicsContext::empty();
        let force = TwoBody::earth();
        let s0 = leo();
        let dt = Second::new(600.0);
        let step = Second::new(6.0); // 100 steps of 6 s = 600 s

        // STM from variational propagator.
        let (_, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();

        // δy0 = 1e-3 km (1 m) in x — linearization error ~ δy0/r ~ 1e-7.
        let delta_y0 = [1e-3_f64, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s_pert = perturb(&s0, &delta_y0);

        // Both nominal and perturbed use the same standard RK4 propagator.
        let y_nom = rk4_propagate(&force, s0, step, 100, &ctx).unwrap();
        let y_pert = rk4_propagate(&force, s_pert, step, 100, &ctx).unwrap();

        let actual_delta: [f64; 6] = {
            let n = state_vec(&y_nom);
            let p = state_vec(&y_pert);
            [
                p[0] - n[0],
                p[1] - n[1],
                p[2] - n[2],
                p[3] - n[3],
                p[4] - n[4],
                p[5] - n[5],
            ]
        };
        let predicted_delta = phi_apply(&phi, &delta_y0);

        let rel_err = max_rel_err(&actual_delta, &predicted_delta);
        assert!(
            rel_err < 1e-6,
            "two-body STM FD check: max relative error = {rel_err:.3e} (tolerance 1e-6)"
        );
    }

    /// J2 STM: same finite-difference validation.
    #[test]
    fn j2_stm_agrees_with_fd() {
        let ctx = DynamicsContext::empty();
        let force = J2::earth();
        let s0 = leo();
        let dt = Second::new(600.0);
        let step = Second::new(6.0);

        let (_, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();

        // δy0 = 1e-3 km (1 m) in y.
        let delta_y0 = [0.0_f64, 1e-3, 0.0, 0.0, 0.0, 0.0];
        let s_pert = perturb(&s0, &delta_y0);

        let y_nom = rk4_propagate(&force, s0, step, 100, &ctx).unwrap();
        let y_pert = rk4_propagate(&force, s_pert, step, 100, &ctx).unwrap();

        let actual_delta: [f64; 6] = {
            let n = state_vec(&y_nom);
            let p = state_vec(&y_pert);
            [
                p[0] - n[0],
                p[1] - n[1],
                p[2] - n[2],
                p[3] - n[3],
                p[4] - n[4],
                p[5] - n[5],
            ]
        };
        let predicted_delta = phi_apply(&phi, &delta_y0);

        let rel_err = max_rel_err(&actual_delta, &predicted_delta);
        assert!(
            rel_err < 1e-6,
            "J2 STM FD check: max relative error = {rel_err:.3e} (tolerance 1e-6)"
        );
    }

    /// Composite (TwoBody + J2) STM agrees with FD within 1e-5 relative.
    #[test]
    fn composite_stm_agrees_with_fd() {
        let ctx = DynamicsContext::empty();
        let force = CompositeForce::empty()
            .push(Box::new(TwoBody::earth()))
            .push(Box::new(J2::earth()));
        let s0 = leo();
        let dt = Second::new(600.0);
        let step = Second::new(6.0);

        let (_, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();

        // δy0 = 1e-3 km (1 m) in z.
        let delta_y0 = [0.0_f64, 0.0, 1e-3, 0.0, 0.0, 0.0];
        let s_pert = perturb(&s0, &delta_y0);

        let y_nom = rk4_propagate(&force, s0, step, 100, &ctx).unwrap();
        let y_pert = rk4_propagate(&force, s_pert, step, 100, &ctx).unwrap();

        let actual_delta: [f64; 6] = {
            let n = state_vec(&y_nom);
            let p = state_vec(&y_pert);
            [
                p[0] - n[0],
                p[1] - n[1],
                p[2] - n[2],
                p[3] - n[3],
                p[4] - n[4],
                p[5] - n[5],
            ]
        };
        let predicted_delta = phi_apply(&phi, &delta_y0);

        let rel_err = max_rel_err(&actual_delta, &predicted_delta);
        assert!(
            rel_err < 1e-5,
            "composite STM FD check: max relative error = {rel_err:.3e} (tolerance 1e-5)"
        );
    }

    /// Φ must have determinant ≈ 1.0 (Liouville / symplectic, conservative forces).
    #[test]
    fn stm_determinant_near_one_two_body() {
        let ctx = DynamicsContext::empty();
        let force = TwoBody::earth();
        let s0 = leo();
        let dt = Second::new(600.0);

        let (_, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();
        let d = det6(*phi.as_array());
        assert!(
            (d - 1.0).abs() < 1e-8,
            "det(Φ) = {d:.10e}, expected ≈ 1.0 (tolerance 1e-8)"
        );
    }

    /// Φ determinant for J2 (also conservative) should also be ≈ 1.0.
    #[test]
    fn stm_determinant_near_one_j2() {
        let ctx = DynamicsContext::empty();
        let force = J2::earth();
        let s0 = leo();
        let dt = Second::new(600.0);

        let (_, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();
        let d = det6(*phi.as_array());
        assert!(
            (d - 1.0).abs() < 1e-8,
            "J2 det(Φ) = {d:.10e}, expected ≈ 1.0 (tolerance 1e-8)"
        );
    }

    /// Determinism: calling propagate_stm twice on the same input gives the same result.
    #[test]
    fn propagate_stm_is_deterministic() {
        let ctx = DynamicsContext::empty();
        let force = TwoBody::earth();
        let s0 = leo();
        let dt = Second::new(300.0);

        let (s1, phi1) = propagate_stm(&force, s0, dt, &ctx).unwrap();
        let (s2, phi2) = propagate_stm(&force, s0, dt, &ctx).unwrap();

        let v1 = state_vec(&s1);
        let v2 = state_vec(&s2);
        for i in 0..6 {
            assert_eq!(v1[i], v2[i], "orbit state component {i} not deterministic");
        }
        let m1 = phi1.as_array();
        let m2 = phi2.as_array();
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(
                    m1[i][j], m2[i][j],
                    "Φ[{i},{j}] not deterministic: {:.16e} vs {:.16e}",
                    m1[i][j], m2[i][j]
                );
            }
        }
    }
}
