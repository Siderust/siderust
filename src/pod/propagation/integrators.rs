// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Uniform [`Integrator`] trait over the upstream RK4 / DOPRI5 / DOP853
//! integrators.
//!
//! ## Why
//!
//! `principia` exposes [`principia::Stepper`] and [`principia::AdaptiveStepper`]
//! as two *different* traits, because RK4 needs a fixed sub-step and the adaptive
//! Hairer-style integrators carry a tolerance configuration. That distinction
//! matters at the algorithm boundary, but most callers — POD batch windows,
//! thrust-arc propagation, validation harnesses — only want a single
//! `propagate(force, state, dt)` interface that they can hand a runtime-
//! selected integrator to.
//!
//! [`Integrator`] is exactly that thin uniform interface. Each impl
//! delegates to the canonical upstream [`principia::rk4_propagate`] /
//! [`principia::dopri5_propagate`] / [`principia::dop853_propagate`] functions
//! and re-raises any upstream [`crate::astro::dynamics::errors::DynamicsError`]
//! through the crate-local [`DynamicsError`].

use crate::astro::dynamics::{DynamicsContext, OrbitState};
use principia::{dop853_propagate, dopri5_propagate, rk4_propagate};
use qtty::{IntegratorTolerances, Second};

use super::error::DynamicsError;
use crate::pod::force::SiderustAccelerationModel;

/// Uniform integrator interface.
///
/// Implementors map a single `propagate(force, state, dt, ctx)` call to the
/// corresponding upstream propagator. `dt` may be negative for backward
/// propagation when the underlying algorithm supports it (DOPRI5 / DOP853).
///
/// # Example
///
/// ```
/// use siderust::pod::propagation::{Integrator, Dop853Integrator};
/// use siderust::pod::force::{OrbitState, Position, Velocity,
///                          DynamicsContext, TwoBody};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::{IntegratorTolerances, Second};
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let integ = Dop853Integrator {
///     tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
/// };
/// let s1 = integ.propagate(&TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0, Second::new(60.0), &DynamicsContext::empty()).unwrap();
/// assert!((s1.epoch - s0.epoch).value() > 0.0);
/// ```
pub trait Integrator {
    /// Propagate `state` by `dt` seconds under `force` and `ctx`.
    fn propagate<FM: SiderustAccelerationModel>(
        &self,
        force: &FM,
        state: OrbitState,
        dt: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState, DynamicsError>;
}

/// Fixed-step classical Runge-Kutta 4th-order integrator.
///
/// Uses [`principia::rk4_propagate`] under the
/// hood. The configured `step` is treated as a magnitude; the wrapper
/// computes `n_steps = ceil(|dt| / step)` so the *effective* sub-step never
/// exceeds `step` in absolute value. Forward propagation is supported (RK4
/// is fixed-step and not naturally signed; pass a positive `dt`).
///
/// # Example
///
/// ```
/// use siderust::pod::propagation::{Integrator, Rk4Integrator};
/// use siderust::pod::force::{OrbitState, Position, Velocity,
///                          DynamicsContext, TwoBody};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::Second;
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let integ = Rk4Integrator { step: Second::new(10.0) };
/// let s1 = integ.propagate(&TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0, Second::new(60.0), &DynamicsContext::empty()).unwrap();
/// assert!((s1.epoch - s0.epoch).value() > 0.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Rk4Integrator {
    /// Maximum sub-step magnitude in seconds. Must be strictly positive.
    pub step: Second,
}

impl Integrator for Rk4Integrator {
    fn propagate<FM: SiderustAccelerationModel>(
        &self,
        force: &FM,
        state: OrbitState,
        dt: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState, DynamicsError> {
        let dt_s = dt.value();
        if dt_s == 0.0 {
            return Ok(state);
        }
        let step = self.step.value().abs();
        if step <= 0.0 || step.is_nan() {
            return Err(DynamicsError::Upstream(
                crate::astro::dynamics::errors::DynamicsError::InvalidStepRequest {
                    reason: "Rk4Integrator.step must be strictly positive",
                },
            ));
        }
        let n_steps = (dt_s.abs() / step).ceil().max(1.0) as usize;
        let h = Second::new(dt_s / n_steps as f64);
        Ok(rk4_propagate(force, state, h, dt, ctx)?)
    }
}

/// Adaptive Dormand-Prince 5(4) integrator (DOPRI5).
///
/// Wraps [`principia::dopri5_propagate`].
/// Suitable for general-purpose LEO/MEO propagation when 5th-order accuracy
/// is enough.
///
/// # Example
///
/// ```
/// use siderust::pod::propagation::{Integrator, Dopri5Integrator};
/// use qtty::IntegratorTolerances;
/// let _ = Dopri5Integrator {
///     tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Dopri5Integrator {
    /// Adaptive-step tolerances forwarded to the upstream propagator.
    pub tolerances: IntegratorTolerances,
}

impl Integrator for Dopri5Integrator {
    fn propagate<FM: SiderustAccelerationModel>(
        &self,
        force: &FM,
        state: OrbitState,
        dt: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState, DynamicsError> {
        Ok(dopri5_propagate(force, state, dt, self.tolerances, ctx)?)
    }
}

/// Adaptive Hairer DOP853 8th-order integrator.
///
/// Wraps [`principia::dop853_propagate`]. This
/// is the high-precision choice for POD batch windows and finite-burn
/// integration where 1e-9 relative tolerance is the working point.
///
/// # Example
///
/// ```
/// use siderust::pod::propagation::{Integrator, Dop853Integrator};
/// use qtty::IntegratorTolerances;
/// let _ = Dop853Integrator {
///     tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Dop853Integrator {
    /// Adaptive-step tolerances forwarded to the upstream propagator.
    pub tolerances: IntegratorTolerances,
}

impl Integrator for Dop853Integrator {
    fn propagate<FM: SiderustAccelerationModel>(
        &self,
        force: &FM,
        state: OrbitState,
        dt: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState, DynamicsError> {
        Ok(dop853_propagate(force, state, dt, self.tolerances, ctx)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::forces::TwoBody;
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
    fn rk4_zero_dt_is_identity() {
        let integ = Rk4Integrator {
            step: Second::new(10.0),
        };
        let s = integ
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s0(),
                Second::new(0.0),
                &DynamicsContext::empty(),
            )
            .unwrap();
        let dt_total = s.epoch - s0().epoch;
        assert!(dt_total.value().abs() < 1e-15);
    }

    #[test]
    fn rk4_invalid_step_rejected() {
        let integ = Rk4Integrator {
            step: Second::new(0.0),
        };
        let r = integ.propagate(
            &TwoBody::new(crate::astro::dynamics::GM_EARTH),
            s0(),
            Second::new(60.0),
            &DynamicsContext::empty(),
        );
        assert!(matches!(
            r,
            Err(DynamicsError::Upstream(
                crate::astro::dynamics::errors::DynamicsError::InvalidStepRequest { .. }
            ))
        ));
    }

    #[test]
    fn dop853_and_dopri5_agree_on_short_arc() {
        let tol = IntegratorTolerances::uniform(1e-12, 1e-9, 1e-12);
        let dt = Second::new(120.0);
        let ctx = DynamicsContext::empty();
        let a = Dop853Integrator { tolerances: tol }
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s0(),
                dt,
                &ctx,
            )
            .unwrap();
        let b = Dopri5Integrator { tolerances: tol }
            .propagate(
                &TwoBody::new(crate::astro::dynamics::GM_EARTH),
                s0(),
                dt,
                &ctx,
            )
            .unwrap();
        let dx = (a.position.x().value() - b.position.x().value()).abs();
        let dy = (a.position.y().value() - b.position.y().value()).abs();
        let dz = (a.position.z().value() - b.position.z().value()).abs();
        assert!(dx + dy + dz < 1e-3, "dx+dy+dz = {}", dx + dy + dz);
    }
}
