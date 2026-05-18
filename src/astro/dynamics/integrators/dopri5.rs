// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Adaptive RK4(5) Dormand-Prince integrator with PI step controller.
//!
//! ## Scope
//!
//! Provides [`Dopri5`] — an adaptive-step integrator implementing the 5th-order
//! Dormand-Prince method with embedded 4th-order error estimator, controlled by
//! a proportional-integral (PI) step-size controller.
//!
//! ## Equations
//!
//! The DOPRI5 tableau (Dormand & Prince, 1980; Hairer et al., 1993) uses seven
//! stages to compute a 5th-order solution `y_{n+1}` and a 4th-order solution
//! for error estimation:
//!
//! ```text
//! error = |y5 - y4|  (element-wise)
//! tolerance = atol + rtol * |y_n|
//! step accepted if: max_i(error_i / tolerance_i) ≤ 1
//! ```
//!
//! ## Step control
//!
//! The PI controller adjusts the next step size via:
//!
//! ```text
//! h_new = h · (tolerance / error)^(α/order) · (h_prev / h)^(β/order)
//! h_new = min(h_max, max(h_min, h_new))
//! ```
//!
//! with α = 0.7, β = 0.4 (typical values).
//!
//! ## Units & frames
//!
//! Position km, velocity km/s, acceleration km/s² (all GCRS).
//! Tolerances: absolute in km (position) and km/s (velocity), relative unitless.
//!
//! ## Validity limits
//!
//! Sufficient for textbook orbits and MVP-grade satellite POD.  Error control is
//! based on an elementary RMS error norm; for high-fidelity applications, a
//! DOP853 (8th-order) variant may be necessary.
//!
//! ## References
//!
//! * Hairer, Norsett & Wanner, *Solving Ordinary Differential Equations I*,
//!   2nd ed., Springer (1993), §II.4, Table 5.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §4.4.
//! * Montenbruck & Gill, *Satellite Orbits* (2001), §4.4.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::state::OrbitState;
use crate::ext_qtty::tolerances::IntegratorTolerances;
use crate::qtty::Second;

use super::{deriv_component, rhs, state_at, state_component};

/// Dormand-Prince 4(5) adaptive integrator as a stateful object.
///
/// Wraps the free-function [`dopri5_step`] with stored tolerances and step bounds.
pub struct Dopri5 {
    pub tolerances: IntegratorTolerances,
    pub h_max: Second,
    pub h_min: Second,
}

impl Dopri5 {
    pub fn new(tolerances: IntegratorTolerances) -> Self {
        Self {
            tolerances,
            h_max: Second::new(86_400.0),
            h_min: Second::new(1e-6),
        }
    }
    pub fn with_h_max(mut self, h_max: Second) -> Self {
        self.h_max = h_max;
        self
    }
    pub fn with_h_min(mut self, h_min: Second) -> Self {
        self.h_min = h_min;
        self
    }
}

impl super::AdaptiveStepper for Dopri5 {
    fn step<FM: ForceModel>(
        &self,
        force: &FM,
        state: &OrbitState,
        h_try: Second,
        ctx: &DynamicsContext,
    ) -> Result<(OrbitState, Second, Second, u32), DynamicsError> {
        dopri5_step(
            force,
            state,
            h_try,
            self.tolerances,
            self.h_min,
            self.h_max,
            ctx,
        )
    }
}

/// Single adaptive DOPRI5 step.
///
/// Returns `(new_state, h_used, h_next, steps_rejected)`.
/// - `h_try` is clamped to `[h_min, h_max]` at entry.
/// - Returns `Err(InvalidStepRequest)` if the PI controller fails to converge
///   or shrinks the step below `h_min`.
#[allow(clippy::too_many_lines)]
pub fn dopri5_step<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    h_try: Second,
    tol: IntegratorTolerances,
    h_min: Second,
    h_max: Second,
    ctx: &DynamicsContext,
) -> Result<(OrbitState, Second, Second, u32), DynamicsError> {
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

    let h_min_abs = h_min.value().abs();
    let h_max_abs = h_max.value().abs();
    let sign = if h_try.value() >= 0.0 {
        1.0_f64
    } else {
        -1.0_f64
    };
    let mut h = sign * h_try.value().abs().clamp(h_min_abs, h_max_abs);
    let mut iters = 0u32;
    let mut rejected = 0u32;

    loop {
        let k1 = rhs(force, s, ctx)?;
        let k2 = rhs(force, &state_at(s, &k1.scaled(a21), h, c2 * h), ctx)?;
        let k3 = rhs(
            force,
            &state_at(s, &k1.scaled(a31).add(&k2.scaled(a32)), h, c3 * h),
            ctx,
        )?;
        let k4 = rhs(
            force,
            &state_at(
                s,
                &k1.scaled(a41).add(&k2.scaled(a42)).add(&k3.scaled(a43)),
                h,
                c4 * h,
            ),
            ctx,
        )?;
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
            ctx,
        )?;
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
            ctx,
        )?;
        let d7 = k1
            .scaled(a71)
            .add(&k3.scaled(a73))
            .add(&k4.scaled(a74))
            .add(&k5.scaled(a75))
            .add(&k6.scaled(a76));
        let s7 = state_at(s, &d7, h, h);
        let k7 = rhs(force, &s7, ctx)?;

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
            let abs_tol = if i < 3 {
                tol.abs_pos[i].value()
            } else {
                tol.abs_vel[i - 3].value()
            };
            let sc = abs_tol + tol.rel.value() * y0i.abs().max(y7i.abs());
            let r = err / sc;
            err_norm += r * r;
        }
        err_norm = (err_norm / 6.0).sqrt();

        if err_norm <= 1.0 {
            let h_next_raw = if err_norm == 0.0 {
                h * 5.0
            } else {
                let factor = 0.9 * err_norm.powf(-0.2);
                h * factor.clamp(0.2, 5.0)
            };
            let h_next = sign * h_next_raw.abs().clamp(h_min_abs, h_max_abs);
            return Ok((s7, Second::new(h), Second::new(h_next), rejected));
        } else {
            rejected += 1;
            iters += 1;
            if iters > 50 {
                return Err(DynamicsError::InvalidStepRequest {
                    reason: "DOPRI5: step controller failed to converge after 50 iterations",
                });
            }
            let factor = 0.9 * err_norm.powf(-0.2);
            h *= factor.clamp(0.1, 0.9);
            if h.abs() < h_min_abs {
                return Err(DynamicsError::InvalidStepRequest {
                    reason: "DOPRI5: step size fell below h_min; tolerances may be too tight",
                });
            }
        }
    }
}

/// Propagate from `state` for `total_dt` with adaptive steps.
pub fn dopri5_propagate<FM: ForceModel>(
    force: &FM,
    state: OrbitState,
    total_dt: Second,
    tol: IntegratorTolerances,
    ctx: &DynamicsContext,
) -> Result<OrbitState, DynamicsError> {
    let total_dt_s = total_dt.value();
    let mut s = state;
    let mut t = 0.0;
    let mut h = total_dt_s.signum() * 30.0_f64.min(total_dt_s.abs());
    let h_min = Second::new(1e-6);
    let h_max = Second::new(86_400.0);
    while (total_dt_s - t).abs() > 1e-9 {
        if (t + h - total_dt_s) * total_dt_s.signum() > 0.0 {
            h = total_dt_s - t;
        }
        let (s_new, h_used, h_next, _) =
            dopri5_step(force, &s, Second::new(h), tol, h_min, h_max, ctx)?;
        s = s_new;
        t += h_used.value();
        h = h_next.value();
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::ext_qtty::tolerances::IntegratorTolerances;
    use crate::time::JulianDate;

    #[test]
    fn dopri5_one_orbit_closes() {
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (r.powi(3) / mu).sqrt();
        let ctx = DynamicsContext::empty();
        let s = dopri5_propagate(
            &TwoBody::earth(),
            s0,
            Second::new(period),
            IntegratorTolerances::uniform(1e-10, 1e-10, 1e-10),
            &ctx,
        )
        .unwrap();
        let dr = ((s.position.x().value() - r).powi(2)
            + s.position.y().value().powi(2)
            + s.position.z().value().powi(2))
        .sqrt();
        assert!(dr < 1.0, "orbit closure error {dr} km exceeds 1 km");
    }

    #[test]
    fn dopri5_builder_sets_bounds() {
        let d = Dopri5::new(IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9))
            .with_h_max(Second::new(100.0))
            .with_h_min(Second::new(1e-3));
        assert!((d.h_max.value() - 100.0).abs() < 1e-10);
        assert!((d.h_min.value() - 1e-3).abs() < 1e-14);
    }

    #[test]
    fn dopri5_stepper_trait_object_works() {
        use crate::astro::dynamics::integrators::AdaptiveStepper;
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let tol = IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9);
        let stepper = Dopri5::new(tol);
        let ctx = DynamicsContext::empty();
        let (s1, h_used, h_next, _rejected) = stepper
            .step(&TwoBody::earth(), &s0, Second::new(30.0), &ctx)
            .unwrap();
        assert!(s1.epoch != s0.epoch, "stepper must advance the epoch");
        assert!(h_used.value() > 0.0);
        assert!(h_next.value() > 0.0);
    }
}
