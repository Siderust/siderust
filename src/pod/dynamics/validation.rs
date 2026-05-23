// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! STM finite-difference validation harness.
//!
//! The variational STM produced by
//! [`principia::propagate_stm`] is the
//! production code path. This module wraps the canonical
//! "perturb-propagate-compare" validation that asserts
//! `Φ · δy₀ ≈ y(t; y₀ + δy₀) − y(t; y₀)` to within a tight tolerance.
//!
//! For Phase 2.3 the maturity plan requires agreement better than `1e-7`
//! (relative). With a 1-m perturbation on a 7000 km orbit, the second-order
//! linearisation error scales as `(|δy₀|/r)² ≈ 2e-14`, so the achievable
//! floor is bounded by RK4 truncation noise rather than by the
//! linearisation, and `< 1e-7` is reached comfortably (typically `< 1e-8`).
//!
//! The function is force-model-agnostic: pass [`siderust::astro::dynamics::forces::TwoBody`],
//! [`siderust::astro::dynamics::forces::J2`], or any [`ForceModel`] whose
//! analytic [`ForceModel::partials`] are wired up.

use principia::{propagate_stm, rk4_propagate};
use crate::astro::dynamics::{DynamicsContext, OrbitState};
use qtty::Second;

use super::error::DynamicsError;
use crate::pod::dynamics::SiderustAccelerationModel;

/// Maximum relative error between the variational-STM-predicted state
/// perturbation and the directly-propagated finite-difference perturbation.
///
/// Protocol:
///
/// 1. Compute `Φ` over `dt` using [`principia::propagate_stm`].
/// 2. Propagate `state` and `state + δy₀` forward by `dt` with RK4 at
///    `step` sub-step.
/// 3. Compare `y(perturb) − y(nominal)` (the *finite-difference* truth) to
///    `Φ · δy₀` element-wise; return the max relative error
///    `|actual − predicted| / max(|actual|, eps)`.
///
/// `delta_y0` must be expressed as `[δrx, δry, δrz, δvx, δvy, δvz]` in
/// km and km/s respectively.
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::validation::max_rel_stm_predict_error;
/// use siderust::pod::dynamics::{DynamicsContext, OrbitState, Position, Velocity, TwoBody};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::Second;
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let err = max_rel_stm_predict_error(
///     &TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0,
///     Second::new(60.0), Second::new(0.6),
///     [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0],
///     &DynamicsContext::empty(),
/// ).unwrap();
/// assert!(err < 1e-7, "STM predict error {err:.3e} exceeds 1e-7");
/// ```
pub fn max_rel_stm_predict_error<FM: SiderustAccelerationModel>(
    force: &FM,
    state: OrbitState,
    dt: Second,
    step: Second,
    delta_y0: [f64; 6],
    ctx: &DynamicsContext,
) -> Result<f64, DynamicsError> {
    let dt_s = dt.value();
    let step_s = step.value().abs().max(1e-12);
    let n_steps = (dt_s.abs() / step_s).ceil().max(1.0) as usize;
    let h = Second::new(dt_s / n_steps as f64);

    let (_, phi) = propagate_stm(force, state, dt, ctx)?;
    let phi_arr = *phi.as_array();

    let nominal = rk4_propagate(force, state, h, dt, ctx)?;
    let perturbed = rk4_propagate(force, perturb_state(&state, &delta_y0), h, dt, ctx)?;

    let n = state_vec(&nominal);
    let p = state_vec(&perturbed);
    let actual = [
        p[0] - n[0],
        p[1] - n[1],
        p[2] - n[2],
        p[3] - n[3],
        p[4] - n[4],
        p[5] - n[5],
    ];
    let predicted = phi_apply(&phi_arr, &delta_y0);

    let mut max = 0.0_f64;
    let scale = actual
        .iter()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max)
        .max(1e-30);
    for i in 0..6 {
        let r = (actual[i] - predicted[i]).abs() / scale;
        if r > max {
            max = r;
        }
    }
    Ok(max)
}

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

fn perturb_state(s: &OrbitState, d: &[f64; 6]) -> OrbitState {
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    OrbitState::new(
        s.epoch,
        Position::<GCRS>::new(
            s.position.x().value() + d[0],
            s.position.y().value() + d[1],
            s.position.z().value() + d[2],
        ),
        Velocity::<GCRS>::new(
            s.velocity.x().value() + d[3],
            s.velocity.y().value() + d[4],
            s.velocity.z().value() + d[5],
        ),
    )
}

fn phi_apply(phi: &[[f64; 6]; 6], dy0: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, row) in phi.iter().enumerate() {
        let mut s = 0.0;
        for (j, v) in row.iter().enumerate() {
            s += v * dy0[j];
        }
        out[i] = s;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::forces::{TwoBody, J2};
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    fn s0() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(7_000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
        )
    }

    #[test]
    fn two_body_variational_matches_finite_diff() {
        let err = max_rel_stm_predict_error(
            &TwoBody::new(crate::astro::dynamics::GM_EARTH),
            s0(),
            Second::new(60.0),
            Second::new(0.6),
            [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0],
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!(
            err < 1e-7,
            "two-body STM predict err {err:.3e} exceeds 1e-7"
        );
    }

    #[test]
    fn j2_variational_matches_finite_diff() {
        let err = max_rel_stm_predict_error(
            &J2::new(crate::astro::dynamics::GM_EARTH, crate::astro::dynamics::R_EARTH, crate::astro::dynamics::EARTH_J2),
            s0(),
            Second::new(60.0),
            Second::new(0.6),
            [0.0, 1e-3, 0.0, 0.0, 0.0, 0.0],
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!(err < 1e-7, "J2 STM predict err {err:.3e} exceeds 1e-7");
    }
}
