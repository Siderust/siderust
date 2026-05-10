// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Adaptive RK4(5) Dormand-Prince integrator with PI step controller.
//!
//! Sufficient for textbook orbits and MVP-grade satellite POD.  A
//! higher-order DOP853 variant is left for follow-up work.

use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::state::{OrbitState, StateDerivative};
use crate::time::JulianDate;

const SEC_PER_DAY: f64 = 86_400.0;

/// Tolerances for adaptive integration.
#[derive(Debug, Clone, Copy)]
pub struct Tolerance {
    /// Relative tolerance.
    pub rel: f64,
    /// Absolute tolerance (km, km/s).
    pub abs: f64,
}

impl Default for Tolerance {
    fn default() -> Self {
        Self { rel: 1e-10, abs: 1e-10 }
    }
}

#[inline]
fn state_component(s: &OrbitState, i: usize) -> f64 {
    match i {
        0 => s.position.x().value(),
        1 => s.position.y().value(),
        2 => s.position.z().value(),
        3 => s.velocity.x().value(),
        4 => s.velocity.y().value(),
        5 => s.velocity.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
fn deriv_component(d: &StateDerivative, i: usize) -> f64 {
    match i {
        0 => d.vel.x().value(),
        1 => d.vel.y().value(),
        2 => d.vel.z().value(),
        3 => d.acc.x().value(),
        4 => d.acc.y().value(),
        5 => d.acc.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
fn rhs<F: ForceModel>(force: &F, s: &OrbitState) -> StateDerivative {
    StateDerivative::new(s.velocity, force.acceleration(s))
}

#[inline]
fn state_at(s: &OrbitState, d: &StateDerivative, h: f64, dt: f64) -> OrbitState {
    let jd = JulianDate::new(s.epoch_tt.jd_value() + dt / SEC_PER_DAY);
    let advanced = s.advance(d, h);
    OrbitState::new(jd, advanced.position, advanced.velocity)
}

/// Single adaptive DOPRI5 step. Returns `(new_state, h_used, h_next)`.
#[allow(clippy::too_many_lines)]
pub fn dopri5_step<F: ForceModel>(
    force: &F,
    s: &OrbitState,
    h_try: f64,
    tol: Tolerance,
) -> (OrbitState, f64, f64) {
    let c2 = 1.0 / 5.0;
    let c3 = 3.0 / 10.0;
    let c4 = 4.0 / 5.0;
    let c5 = 8.0 / 9.0;

    let a21 = 1.0 / 5.0;
    let a31 = 3.0 / 40.0;
    let a32 = 9.0 / 40.0;
    let a41 = 44.0 / 45.0;
    let a42 = -56.0 / 15.0;
    let a43 = 32.0 / 9.0;
    let a51 = 19_372.0 / 6_561.0;
    let a52 = -25_360.0 / 2_187.0;
    let a53 = 64_448.0 / 6_561.0;
    let a54 = -212.0 / 729.0;
    let a61 = 9_017.0 / 3_168.0;
    let a62 = -355.0 / 33.0;
    let a63 = 46_732.0 / 5_247.0;
    let a64 = 49.0 / 176.0;
    let a65 = -5_103.0 / 18_656.0;
    let a71 = 35.0 / 384.0;
    let a73 = 500.0 / 1_113.0;
    let a74 = 125.0 / 192.0;
    let a75 = -2_187.0 / 6_784.0;
    let a76 = 11.0 / 84.0;

    let e1 = 71.0 / 57_600.0;
    let e3 = -71.0 / 16_695.0;
    let e4 = 71.0 / 1_920.0;
    let e5 = -17_253.0 / 339_200.0;
    let e6 = 22.0 / 525.0;
    let e7 = -1.0 / 40.0;

    let mut h = h_try;

    loop {
        let k1 = rhs(force, s);
        let k2 = rhs(force, &state_at(s, &k1.scaled(a21), h, c2 * h));
        let k3 = rhs(force, &state_at(s, &k1.scaled(a31).add(&k2.scaled(a32)), h, c3 * h));
        let k4 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a41).add(&k2.scaled(a42)).add(&k3.scaled(a43)),
                h,
                c4 * h,
            ),
        );
        let k5 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a51)
                    .add(&k2.scaled(a52))
                    .add(&k3.scaled(a53))
                    .add(&k4.scaled(a54)),
                h,
                c5 * h,
            ),
        );
        let k6 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a61)
                    .add(&k2.scaled(a62))
                    .add(&k3.scaled(a63))
                    .add(&k4.scaled(a64))
                    .add(&k5.scaled(a65)),
                h,
                h,
            ),
        );
        let d7 = k1
            .scaled(a71)
            .add(&k3.scaled(a73))
            .add(&k4.scaled(a74))
            .add(&k5.scaled(a75))
            .add(&k6.scaled(a76));
        let s7 = state_at(s, &d7, h, h);
        let k7 = rhs(force, &s7);

        let err_d = k1
            .scaled(e1)
            .add(&k3.scaled(e3))
            .add(&k4.scaled(e4))
            .add(&k5.scaled(e5))
            .add(&k6.scaled(e6))
            .add(&k7.scaled(e7));
        let mut err_norm = 0.0;
        for i in 0..6 {
            let err = h * deriv_component(&err_d, i);
            let y0i = state_component(s, i);
            let y7i = state_component(&s7, i);
            let sc = tol.abs + tol.rel * y0i.abs().max(y7i.abs());
            let r = err / sc;
            err_norm += r * r;
        }
        err_norm = (err_norm / 6.0).sqrt();

        if err_norm <= 1.0 || h.abs() < 1e-9 {
            let h_next = if err_norm == 0.0 {
                h * 5.0
            } else {
                let factor = 0.9 * err_norm.powf(-0.2);
                h * factor.clamp(0.2, 5.0)
            };
            return (s7, h, h_next);
        } else {
            let factor = 0.9 * err_norm.powf(-0.2);
            h *= factor.clamp(0.1, 0.9);
        }
    }
}

/// Propagate from `state` for `total_dt_s` seconds with adaptive steps.
pub fn dopri5_propagate<F: ForceModel>(
    force: &F,
    state: OrbitState,
    total_dt_s: f64,
    tol: Tolerance,
) -> OrbitState {
    let mut s = state;
    let mut t = 0.0;
    let mut h = total_dt_s.signum() * 30.0_f64.min(total_dt_s.abs());
    while (total_dt_s - t).abs() > 1e-9 {
        if (t + h - total_dt_s) * total_dt_s.signum() > 0.0 {
            h = total_dt_s - t;
        }
        let (s_new, h_used, h_next) = dopri5_step(force, &s, h, tol);
        s = s_new;
        t += h_used;
        h = h_next;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;

    #[test]
    fn dopri5_one_orbit_closes() {
        let mu: f64 = 398_600.4418;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (r.powi(3) / mu).sqrt();
        let s = dopri5_propagate(&TwoBody::earth(), s0, period, Tolerance::default());
        let dr = ((s.position.x().value() - r).powi(2)
            + s.position.y().value().powi(2)
            + s.position.z().value().powi(2))
        .sqrt();
        assert!(dr < 1.0, "orbit closure error {dr} km exceeds 1 km");
    }
}
