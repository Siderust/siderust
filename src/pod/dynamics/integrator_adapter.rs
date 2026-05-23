// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed integrator adapter.
//!
//! Re-exports of [`siderust::pod::dynamics::Integrator`] (and the concrete
//! Rk4/Dopri5/Dop853 wrappers) live in the crate root. This module provides
//! the POD-level convenience helpers that take a [`SpacecraftState`] (or
//! [`OrbitState`]) plus a registry-built composite force model and return
//! the propagated state.
//!
//! No new integrator math is invented here — every call delegates to the
//! upstream propagator.

use crate::astro::dynamics::{DynamicsContext, OrbitState, SpacecraftState};
use qtty::Second;

use crate::pod::dynamics::{Integrator, SiderustAccelerationModel};

use super::pod_error::PodDynamicsError;

/// Propagate a bare [`OrbitState`] forward by `dt` under `force` using the
/// supplied [`Integrator`].
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::{propagate_orbit, Rk4Integrator};
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::{OrbitState, DynamicsContext};
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::Second;
///
/// let s0 = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let s1 = propagate_orbit(
///     &Rk4Integrator { step: Second::new(10.0) },
///     &TwoBody::new(siderust::astro::dynamics::GM_EARTH), s0, Second::new(60.0), &DynamicsContext::empty(),
/// ).unwrap();
/// assert!((s1.epoch - s0.epoch).value() > 0.0);
/// ```
pub fn propagate_orbit<I, F>(
    integrator: &I,
    force: &F,
    state: OrbitState,
    dt: Second,
    ctx: &DynamicsContext,
) -> Result<OrbitState, PodDynamicsError>
where
    I: Integrator,
    F: SiderustAccelerationModel,
{
    Ok(integrator.propagate(force, state, dt, ctx)?)
}

/// Propagate a [`SpacecraftState`] (orbit + properties) forward by `dt`.
///
/// Spacecraft properties (`C_D`, `C_R`, `A/m`, mass) are passed through
/// unchanged — finite-burn / mass-flow integration is the responsibility of
/// the upstream `siderust::pod::dynamics::thrust` module, not this adapter.
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::{propagate_spacecraft, Rk4Integrator};
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::{
///     DynamicsContext, OrbitState, SpacecraftProperties, SpacecraftState,
/// };
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
/// use qtty::Second;
///
/// let orbit = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
/// );
/// let sc = SpacecraftState { orbit, properties: SpacecraftProperties::demo_leo() };
/// let sc1 = propagate_spacecraft(
///     &Rk4Integrator { step: Second::new(10.0) },
///     &TwoBody::new(siderust::astro::dynamics::GM_EARTH), sc, Second::new(60.0), &DynamicsContext::empty(),
/// ).unwrap();
/// assert_eq!(sc1.properties, sc.properties);
/// ```
pub fn propagate_spacecraft<I, F>(
    integrator: &I,
    force: &F,
    state: SpacecraftState,
    dt: Second,
    ctx: &DynamicsContext,
) -> Result<SpacecraftState, PodDynamicsError>
where
    I: Integrator,
    F: SiderustAccelerationModel,
{
    let orbit = integrator.propagate(force, state.orbit, dt, ctx)?;
    Ok(SpacecraftState {
        orbit,
        properties: state.properties,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pod::dynamics::Rk4Integrator;
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
    fn orbit_advances() {
        let integ = Rk4Integrator {
            step: Second::new(10.0),
        };
        let s1 = propagate_orbit(
            &integ,
            &TwoBody::new(crate::astro::dynamics::GM_EARTH),
            s0(),
            Second::new(60.0),
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!((s1.epoch - s0().epoch).value() > 0.0);
    }
}
