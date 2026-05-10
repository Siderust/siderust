// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Numerical state-transition matrices for Cartesian orbital states.
//!
//! Computes the Jacobian Φ(t,t₀) = ∂x(t)/∂x(t₀) by central finite differences:
//! each of the six state components is perturbed by ±h (relative, `h = 1e-6 · max(|x₀ⱼ|, 1)`),
//! propagated with the supplied force model and RK4 integrator, and the column is formed
//! from the central-difference quotient.
//!
//! This avoids coding analytic variational equations per force model and is
//! accurate enough for batch POD windows up to a few hours.

use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::integrators::{rk4_propagate, rk4_propagate_series};
use crate::astro::dynamics::{OrbitState, Position, Velocity};
use crate::coordinates::frames::GCRS;

/// A 6×6 row-major state-transition (or Jacobian) matrix.
pub type Stm6 = [[f64; 6]; 6];

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
        s.epoch_tt,
        Position::<GCRS>::new(rx, ry, rz),
        Velocity::<GCRS>::new(vx, vy, vz),
    )
}

/// Finite-difference state-transition matrix Φ(t, t₀).
///
/// Returns a [`Stm6`] row-major Jacobian `∂x(t)/∂x(t₀)`.
/// `dt_s` is the RK4 step size in seconds; `n_steps` is the total number of steps.
pub fn finite_diff_stm<F: ForceModel>(
    force: &F,
    s0: OrbitState,
    dt_s: f64,
    n_steps: usize,
) -> Stm6 {
    let mut stm = [[0.0; 6]; 6];
    for j in 0..6 {
        let x0j = state_component(&s0, j);
        let h = 1e-6 * x0j.abs().max(1.0);
        let s_plus = rk4_propagate(force, perturb_component(&s0, j, h), dt_s, n_steps);
        let s_minus = rk4_propagate(force, perturb_component(&s0, j, -h), dt_s, n_steps);
        for i in 0..6 {
            stm[i][j] =
                (state_component(&s_plus, i) - state_component(&s_minus, i)) / (2.0 * h);
        }
    }
    stm
}

/// Compute Φ(tₖ, t₀) at every step `k = 0..=n_steps`.
///
/// Φ₀ is the identity. The 12 perturbation propagations are shared across all
/// output epochs, which is far cheaper than calling [`finite_diff_stm`] per epoch.
pub fn finite_diff_stm_series<F: ForceModel>(
    force: &F,
    s0: OrbitState,
    dt_s: f64,
    n_steps: usize,
) -> Vec<Stm6> {
    let mut perturbed: [[Vec<[f64; 6]>; 2]; 6] = Default::default();
    let mut hs = [0.0_f64; 6];
    for j in 0..6 {
        let x0j = state_component(&s0, j);
        let h = 1e-6 * x0j.abs().max(1.0);
        hs[j] = h;
        for (sign_idx, sign) in [-1.0_f64, 1.0_f64].iter().enumerate() {
            let sp = perturb_component(&s0, j, sign * h);
            perturbed[j][sign_idx] = rk4_propagate_series(force, sp, dt_s, n_steps)
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
        let mut phi = [[0.0; 6]; 6];
        for j in 0..6 {
            let plus = perturbed[j][1][k];
            let minus = perturbed[j][0][k];
            for i in 0..6 {
                phi[i][j] = (plus[i] - minus[i]) / (2.0 * hs[j]);
            }
        }
        out.push(phi);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::forces::{ForceModel, TwoBody};
    use crate::time::JulianDate;

    #[test]
    fn two_body_stm_is_identity_at_zero_steps() {
        let s = OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let phi = finite_diff_stm(&TwoBody::earth(), s, 1.0, 0);
        for i in 0..6 {
            for j in 0..6 {
                let target = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (phi[i][j] - target).abs() < 1e-9,
                    "phi[{i}][{j}] = {} (expected {target})",
                    phi[i][j]
                );
            }
        }
    }

    #[test]
    fn stm_series_first_element_is_identity() {
        let s = OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let series = finite_diff_stm_series(&TwoBody::earth(), s, 1.0, 2);
        assert_eq!(series.len(), 3);
        let phi0 = series[0];
        for i in 0..6 {
            for j in 0..6 {
                let target = if i == j { 1.0 } else { 0.0 };
                assert!((phi0[i][j] - target).abs() < 1e-9);
            }
        }
    }
}
