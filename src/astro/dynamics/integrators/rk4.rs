//! Classical Runge-Kutta 4 fixed-step integrator for [`OrbitState`].
//!
//! ## Scope
//!
//! Provides [`Rk4`] — a zero-sized fixed-step integrator implementing the
//! classical 4th-order Runge-Kutta scheme (RK4).
//!
//! ## Equations
//!
//! The RK4 step from state `y_n` over interval `dt` is:
//!
//! ```text
//! k1 = f(y_n)
//! k2 = f(y_n + dt/2 · k1)
//! k3 = f(y_n + dt/2 · k2)
//! k4 = f(y_n + dt · k3)
//! y_{n+1} = y_n + dt/6 · (k1 + 2k2 + 2k3 + k4)
//! ```
//!
//! where `f(y) = [v; a]` is the state derivative.
//!
//! ## Units & frames
//!
//! Position km, velocity km/s, acceleration km/s² (all GCRS).
//! Time step in seconds.
//!
//! ## Validity limits
//!
//! RK4 is 4th-order accurate per step: `O(dt⁵)` local truncation error.
//! Step size must be chosen small enough to resolve force-model dynamics
//! (typically 30–60 seconds for LEO with J2+drag).  No adaptive step control.
//!
//! ## References
//!
//! * Hairer, Norsett & Wanner, *Solving Ordinary Differential Equations I*,
//!   §II.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §4.4.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::state::{OrbitState, StateDerivative};
use crate::qtty::Second;

/// Classical RK4 as a stateful zero-sized object.
pub struct Rk4;

impl super::FixedStepper for Rk4 {
    fn step<FM: ForceModel>(
        &self,
        force: &FM,
        state: &OrbitState,
        h: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState, DynamicsError> {
        rk4_step(force, state, h, ctx)
    }
}

/// Step the orbit state by `dt` using RK4 with the given force model.
///
/// The returned state has its epoch advanced by `dt`.
pub fn rk4_step<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    dt: Second,
    ctx: &DynamicsContext,
) -> Result<OrbitState, DynamicsError> {
    let dt_s = dt.value();
    let half = Second::new(dt_s * 0.5);
    let k1 = derivative(force, s, ctx)?;
    let s1 = s.advance(&k1, half);
    let k2 = derivative(force, &s1, ctx)?;
    let s2 = s.advance(&k2, half);
    let k3 = derivative(force, &s2, ctx)?;
    let s3 = s.advance(&k3, dt);
    let k4 = derivative(force, &s3, ctx)?;

    let combined = StateDerivative::rk4_combine(&k1, &k2, &k3, &k4);
    let new_epoch = s.epoch + dt;
    let advanced = s.advance(&combined, dt);
    Ok(OrbitState {
        epoch: new_epoch,
        position: advanced.position,
        velocity: advanced.velocity,
    })
}

#[inline]
fn derivative<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    ctx: &DynamicsContext,
) -> Result<StateDerivative, DynamicsError> {
    Ok(StateDerivative::new(
        s.velocity,
        force.acceleration(s, ctx)?,
    ))
}

/// Propagate `state` over `n_steps` of `dt` each.
pub fn rk4_propagate<FM: ForceModel>(
    force: &FM,
    state: OrbitState,
    dt: Second,
    n_steps: usize,
    ctx: &DynamicsContext,
) -> Result<OrbitState, DynamicsError> {
    let mut s = state;
    for _ in 0..n_steps {
        s = rk4_step(force, &s, dt, ctx)?;
    }
    Ok(s)
}

/// Propagate and collect the state at the start and after each step.
/// Output length is `n_steps + 1`.
pub fn rk4_propagate_series<FM: ForceModel>(
    force: &FM,
    state: OrbitState,
    dt: Second,
    n_steps: usize,
    ctx: &DynamicsContext,
) -> Result<Vec<OrbitState>, DynamicsError> {
    let mut out = Vec::with_capacity(n_steps + 1);
    out.push(state);
    let mut s = state;
    for _ in 0..n_steps {
        s = rk4_step(force, &s, dt, ctx)?;
        out.push(s);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    #[test]
    fn rk4_one_orbit_closes() {
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (r.powi(3) / mu).sqrt();
        let n_steps = 4_000;
        let dt = Second::new(period / n_steps as f64);
        let ctx = DynamicsContext::empty();
        let s = rk4_propagate(&TwoBody::earth(), s0, dt, n_steps, &ctx).unwrap();
        let dr = ((s.position.x().value() - r).powi(2)
            + s.position.y().value().powi(2)
            + s.position.z().value().powi(2))
        .sqrt();
        assert!(dr < 5.0, "orbit closure error {dr} km exceeds 5 km");
    }

    #[test]
    fn rk4_propagate_series_length_is_n_steps_plus_one() {
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let ctx = DynamicsContext::empty();
        let series =
            rk4_propagate_series(&TwoBody::earth(), s0, Second::new(60.0), 5, &ctx).unwrap();
        assert_eq!(series.len(), 6, "series must have n_steps+1 elements");
        assert_eq!(series[0], s0, "first element must be the initial state");
    }

    #[test]
    fn rk4_fixed_stepper_advances_epoch() {
        use super::super::FixedStepper;
        let mu: f64 = 398_600.441_8;
        let r: f64 = 7_000.0;
        let v: f64 = (mu / r).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let ctx = DynamicsContext::empty();
        let s1 = Rk4
            .step(&TwoBody::earth(), &s0, Second::new(60.0), &ctx)
            .unwrap();
        assert!(s1.epoch != s0.epoch, "Rk4::step must advance the epoch");
    }
}
