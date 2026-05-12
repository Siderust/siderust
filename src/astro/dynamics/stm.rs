// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon
#![allow(deprecated)]

//! Numerical state-transition matrices for Cartesian orbital states.
//!
//! Computes the Jacobian Φ(t,t₀) = ∂x(t)/∂x(t₀) by central finite differences:
//! each of the six state components is perturbed by ±h (relative, `h = 1e-6 · max(|x₀ⱼ|, 1)`),
//! propagated with the supplied force model and RK4 integrator, and the column is formed
//! from the central-difference quotient.
//!
//! This avoids coding analytic variational equations per force model and is
//! accurate enough for batch POD windows up to a few hours.
//!
//! The result is wrapped in [`StateTransition<F>`], which stores the 6×6 Jacobian as four
//! typed 3×3 blocks:
//!
//! ```text
//! Φ = [ dr_dr  dr_dv ]
//!     [ dv_dr  dv_dv ]
//! ```
//!
//! Call [`StateTransition::to_row_major`] for a dense `[[f64; 6]; 6]` export.

use affn::matrix3::FrameMatrix3;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::integrators::{rk4_propagate, rk4_propagate_series};
use crate::astro::dynamics::{OrbitState, Position, Velocity};
use crate::coordinates::frames::GCRS;
use crate::qtty::Second;

/// Frame-tagged 6×6 state-transition matrix stored as four 3×3 blocks.
///
/// # Deprecation
///
/// This type is the finite-difference STM, preserved for validation only.
/// New code should use [`siderust::astro::dynamics::variational::StateTransitionMatrix`],
/// which is the variational (analytic) STM.
#[deprecated(
    note = "Use `astro::dynamics::variational::StateTransitionMatrix` (`FrameMatrix6<F>`) \
            for production code; this type is preserved for validation only."
)]
#[derive(Debug, Clone, Copy)]
pub struct StateTransition<F> {
    dr_dr: FrameMatrix3<F>,
    dr_dv: FrameMatrix3<F>,
    dv_dr: FrameMatrix3<F>,
    dv_dv: FrameMatrix3<F>,
}

impl<F> StateTransition<F> {
    /// Construct from four 3×3 blocks.
    #[inline]
    pub fn from_blocks(
        dr_dr: FrameMatrix3<F>,
        dr_dv: FrameMatrix3<F>,
        dv_dr: FrameMatrix3<F>,
        dv_dv: FrameMatrix3<F>,
    ) -> Self {
        Self {
            dr_dr,
            dr_dv,
            dv_dr,
            dv_dv,
        }
    }

    /// Identity state-transition: `dr_dr = I`, `dv_dv = I`, off-diagonal = 0.
    pub fn identity() -> Self {
        let eye = || FrameMatrix3::from_array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        Self {
            dr_dr: eye(),
            dr_dv: FrameMatrix3::zero(),
            dv_dr: FrameMatrix3::zero(),
            dv_dv: eye(),
        }
    }

    /// Position–position Jacobian block (`∂r(t)/∂r(t₀)`, dimensionless).
    #[inline]
    pub fn dr_dr(&self) -> &FrameMatrix3<F> {
        &self.dr_dr
    }

    /// Position–velocity Jacobian block (`∂r(t)/∂v(t₀)`, units: time).
    #[inline]
    pub fn dr_dv(&self) -> &FrameMatrix3<F> {
        &self.dr_dv
    }

    /// Velocity–position Jacobian block (`∂v(t)/∂r(t₀)`, units: 1/time).
    #[inline]
    pub fn dv_dr(&self) -> &FrameMatrix3<F> {
        &self.dv_dr
    }

    /// Velocity–velocity Jacobian block (`∂v(t)/∂v(t₀)`, dimensionless).
    #[inline]
    pub fn dv_dv(&self) -> &FrameMatrix3<F> {
        &self.dv_dv
    }

    /// Assemble the full 6×6 row-major matrix for serialisation or numerical routines.
    ///
    /// Layout: rows/columns `0..2` are position, `3..5` are velocity.
    pub fn to_row_major(&self) -> [[f64; 6]; 6] {
        let dr_dr = self.dr_dr.as_array();
        let dr_dv = self.dr_dv.as_array();
        let dv_dr = self.dv_dr.as_array();
        let dv_dv = self.dv_dv.as_array();
        let mut out = [[0.0_f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = dr_dr[i][j];
                out[i][j + 3] = dr_dv[i][j];
                out[i + 3][j] = dv_dr[i][j];
                out[i + 3][j + 3] = dv_dv[i][j];
            }
        }
        out
    }

    /// Re-tag as a different frame without changing data.
    #[inline]
    pub fn relabel<G>(self) -> StateTransition<G> {
        StateTransition {
            dr_dr: self.dr_dr.relabel::<G>(),
            dr_dv: self.dr_dv.relabel::<G>(),
            dv_dr: self.dv_dr.relabel::<G>(),
            dv_dv: self.dv_dv.relabel::<G>(),
        }
    }
}

// =============================================================================
// Private helpers (shared with finite_diff_stm functions)
// =============================================================================

fn state_component(s: &OrbitState, j: usize) -> f64 {
    match j {
        0 => s.position.x().value(),
        1 => s.position.y().value(),
        2 => s.position.z().value(),
        3 => s.velocity.x().value(),
        4 => s.velocity.y().value(),
        5 => s.velocity.z().value(),
        _ => panic!("index out of range"),
    }
}

fn perturb_component(s: &OrbitState, j: usize, delta: f64) -> OrbitState {
    let rx = s.position.x().value() + if j == 0 { delta } else { 0.0 };
    let ry = s.position.y().value() + if j == 1 { delta } else { 0.0 };
    let rz = s.position.z().value() + if j == 2 { delta } else { 0.0 };
    let vx = s.velocity.x().value() + if j == 3 { delta } else { 0.0 };
    let vy = s.velocity.y().value() + if j == 4 { delta } else { 0.0 };
    let vz = s.velocity.z().value() + if j == 5 { delta } else { 0.0 };
    OrbitState::new(
        s.epoch,
        Position::<GCRS>::new(rx, ry, rz),
        Velocity::<GCRS>::new(vx, vy, vz),
    )
}

fn fill_stm_from_raw(raw: &[[f64; 6]; 6]) -> StateTransition<GCRS> {
    let mut dr_dr = [[0.0; 3]; 3];
    let mut dr_dv = [[0.0; 3]; 3];
    let mut dv_dr = [[0.0; 3]; 3];
    let mut dv_dv = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            dr_dr[i][j] = raw[i][j];
            dr_dv[i][j] = raw[i][j + 3];
            dv_dr[i][j] = raw[i + 3][j];
            dv_dv[i][j] = raw[i + 3][j + 3];
        }
    }
    StateTransition::from_blocks(
        FrameMatrix3::from_array(dr_dr),
        FrameMatrix3::from_array(dr_dv),
        FrameMatrix3::from_array(dv_dr),
        FrameMatrix3::from_array(dv_dv),
    )
}

// =============================================================================
// Public API
// =============================================================================

/// Finite-difference state-transition matrix Φ(t, t₀).
///
/// Returns a [`StateTransition<GCRS>`] with four typed 3×3 blocks.
/// `dt` is the RK4 step size; `n_steps` is the total number of steps.
///
/// # Deprecation
///
/// The finite-difference approach is preserved for **validation only**.
/// New code should use the variational-equations propagator:
///
/// ```rust
/// use siderust::astro::dynamics::variational::propagate_stm;
/// ```
#[deprecated(
    note = "Use `astro::dynamics::variational::propagate_stm` for production code; \
            finite-difference STM is preserved for validation only."
)]
pub fn finite_diff_stm<F: ForceModel>(
    force: &F,
    s0: OrbitState,
    dt: Second,
    n_steps: usize,
    ctx: &DynamicsContext,
) -> Result<StateTransition<GCRS>, DynamicsError> {
    let mut raw = [[0.0; 6]; 6];
    for j in 0..6 {
        let x0j = state_component(&s0, j);
        let h = 1e-6 * x0j.abs().max(1.0);
        let s_plus = rk4_propagate(force, perturb_component(&s0, j, h), dt, n_steps, ctx)?;
        let s_minus = rk4_propagate(force, perturb_component(&s0, j, -h), dt, n_steps, ctx)?;
        for (i, row) in raw.iter_mut().enumerate() {
            row[j] = (state_component(&s_plus, i) - state_component(&s_minus, i)) / (2.0 * h);
        }
    }
    Ok(fill_stm_from_raw(&raw))
}

/// Compute Φ(tₖ, t₀) at every step `k = 0..=n_steps`.
///
/// Φ₀ is the identity. The 12 perturbation propagations are shared across all
/// output epochs, which is far cheaper than calling [`finite_diff_stm`] per epoch.
///
/// # Deprecation
///
/// Prefer [`siderust::astro::dynamics::variational::propagate_stm`] for new code.
/// This series variant has no direct variational equivalent and is retained for
/// validation workflows that need the STM at every intermediate step.
#[deprecated(
    note = "Use `astro::dynamics::variational::propagate_stm` for production code; \
            finite-difference STM is preserved for validation only."
)]
pub fn finite_diff_stm_series<F: ForceModel>(
    force: &F,
    s0: OrbitState,
    dt: Second,
    n_steps: usize,
    ctx: &DynamicsContext,
) -> Result<Vec<StateTransition<GCRS>>, DynamicsError> {
    let mut perturbed: [[Vec<[f64; 6]>; 2]; 6] = Default::default();
    let mut hs = [0.0_f64; 6];
    for j in 0..6 {
        let x0j = state_component(&s0, j);
        let h = 1e-6 * x0j.abs().max(1.0);
        hs[j] = h;
        for (sign_idx, sign) in [-1.0_f64, 1.0_f64].iter().enumerate() {
            let sp = perturb_component(&s0, j, sign * h);
            perturbed[j][sign_idx] = rk4_propagate_series(force, sp, dt, n_steps, ctx)?
                .into_iter()
                .map(|s| {
                    [
                        state_component(&s, 0),
                        state_component(&s, 1),
                        state_component(&s, 2),
                        state_component(&s, 3),
                        state_component(&s, 4),
                        state_component(&s, 5),
                    ]
                })
                .collect();
        }
    }
    let mut out = Vec::with_capacity(n_steps + 1);
    for k in 0..=n_steps {
        let mut raw = [[0.0; 6]; 6];
        for (j, perturb_pair) in perturbed.iter().enumerate() {
            let plus = perturb_pair[1][k];
            let minus = perturb_pair[0][k];
            for (i, row) in raw.iter_mut().enumerate() {
                row[j] = (plus[i] - minus[i]) / (2.0 * hs[j]);
            }
        }
        out.push(fill_stm_from_raw(&raw));
    }
    Ok(out)
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::time::JulianDate;

    fn sample_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn check_identity(phi: &[[f64; 6]; 6], tol: f64) {
        for (i, row) in phi.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let target = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (*value - target).abs() < tol,
                    "phi[{i}][{j}] = {} (expected {target})",
                    value
                );
            }
        }
    }

    #[test]
    fn two_body_stm_is_identity_at_zero_steps() {
        let ctx = DynamicsContext::empty();
        let s = sample_state();
        let phi = finite_diff_stm(&TwoBody::earth(), s, Second::new(1.0), 0, &ctx).unwrap();
        check_identity(&phi.to_row_major(), 1e-9);
    }

    #[test]
    fn identity_to_row_major() {
        let phi = StateTransition::<GCRS>::identity().to_row_major();
        check_identity(&phi, 1e-15);
    }

    #[test]
    fn stm_series_first_element_is_identity() {
        let ctx = DynamicsContext::empty();
        let s = sample_state();
        let series =
            finite_diff_stm_series(&TwoBody::earth(), s, Second::new(1.0), 2, &ctx).unwrap();
        assert_eq!(series.len(), 3);
        check_identity(&series[0].to_row_major(), 1e-9);
    }

    #[test]
    fn from_blocks_assembles_six_by_six() {
        let m = FrameMatrix3::<GCRS>::from_array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        let z = FrameMatrix3::<GCRS>::zero();
        let phi = StateTransition::from_blocks(m, z, z, m);
        let raw = phi.to_row_major();
        assert!((raw[0][0] - 1.0).abs() < 1e-15);
        assert!((raw[0][2] - 3.0).abs() < 1e-15);
        assert!((raw[1][1] - 5.0).abs() < 1e-15);
        assert!((raw[3][3] - 1.0).abs() < 1e-15);
        assert!((raw[5][5] - 9.0).abs() < 1e-15);
        assert!((raw[0][3] - 0.0).abs() < 1e-15);
    }

    #[test]
    fn block_accessors_return_correct_blocks() {
        let phi = StateTransition::<GCRS>::identity();
        let dr_dr = phi.dr_dr().as_array();
        let dr_dv = phi.dr_dv().as_array();
        let dv_dr = phi.dv_dr().as_array();
        let dv_dv = phi.dv_dv().as_array();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dr_dr[i][j] - expected).abs() < 1e-15);
                assert!((dv_dv[i][j] - expected).abs() < 1e-15);
                assert!((dr_dv[i][j]).abs() < 1e-15);
                assert!((dv_dr[i][j]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn relabel_preserves_data() {
        use crate::coordinates::frames::ICRS;
        let phi = StateTransition::<GCRS>::identity();
        let relabeled: StateTransition<ICRS> = phi.relabel::<ICRS>();
        check_identity(&relabeled.to_row_major(), 1e-15);
    }

    #[test]
    fn finite_diff_stm_non_identity_after_propagation() {
        let ctx = DynamicsContext::empty();
        let s = sample_state();
        let phi = finite_diff_stm(&TwoBody::earth(), s, Second::new(60.0), 10, &ctx).unwrap();
        let raw = phi.to_row_major();
        // After propagation, the STM must differ from the identity matrix in the
        // dr/dr off-diagonal (orbital coupling).
        let off_diag_max = raw[0][1].abs() + raw[1][0].abs() + raw[0][3].abs() + raw[1][4].abs();
        assert!(off_diag_max > 1e-6, "expected non-trivial STM, got off-diag mass {off_diag_max}");
    }
}
