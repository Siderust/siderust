// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astrodynamics and POD-dynamics composition layer.
//!
//! ## Contents
//!
//! - [`integrators`] — uniform [`Integrator`] trait over RK4/DOPRI5/DOP853
//! - [`thrust`] — finite-burn thrust-arc physical model
//! - [`low_thrust`] — Tsiolkovsky ΔV bookkeeping
//! - [`validation`] — STM finite-difference validation harness
//! - [`empirical_periodic`] — 1-CPR / 2-CPR periodic empirical accelerations
//! - [`force_config`] — configurable force-model spec
//! - [`integrator_adapter`] — typed integrator adapters
//! - [`process_noise`] — Q-matrix construction for EKF
//! - [`registry`] — named force-model registry
//! - [`thrust_arc`] — thrust-arc parameter declaration
//! - [`variational`] — variational/STM propagator
//!
//! ## Force models
//!
//! Force models implement the canonical
//! [`siderust::astro::dynamics::forces::AccelerationModel`] trait and are
//! composed via [`SiderustCompositeModel`]. The built-in models
//! (`TwoBody`, `J2`, `DragForce`, `CannonballSrp`, …) are re-exported
//! from [`siderust::astro::dynamics`] for convenience.
//! The [`registry`] module provides a string-keyed factory for assembling
//! composites from configuration.

pub mod empirical_periodic;
pub mod error;
pub mod force_config;
pub mod integrator_adapter;
pub mod integrators;
pub mod low_thrust;
pub mod pod_error;
pub mod process_noise;
pub mod registry;
pub mod thrust;
pub mod thrust_arc;
pub mod validation;
pub mod variational;

// Re-exports from the low-level layer (former siderust-dynamics public surface)
pub use error::DynamicsError;
pub use integrators::{Dop853Integrator, Dopri5Integrator, Integrator, Rk4Integrator};
pub use low_thrust::{LowThrustLog, LowThrustRecord};
pub use thrust::{mass_flow_rate, thrust_acceleration, ManeuverError, ThrustArc, G0_M_PER_S2};

// Re-exports from the POD dynamics layer
pub use empirical_periodic::{EmpiricalPeriodicAcceleration, PeriodicHarmonic};
pub use force_config::ForceModelConfig;
pub use integrator_adapter::{propagate_orbit, propagate_spacecraft};
pub use pod_error::PodDynamicsError;
pub use process_noise::{
    GaussMarkovParams, PiecewiseSegment, ProcessNoise, ProcessNoiseModel, WhiteAccelPsd,
};
pub use registry::{ForceModelFactory, ForceModelParams, ForceModelRegistry, ForceModelSpec};
pub use thrust_arc::ThrustArcConfig;
pub use variational::{
    param_partials_central_diff, ParamColumn, ParamStmReport, PropagatedArc, VariationalPropagator,
};

/// Geocentric GCRS Cartesian orbit state (position, velocity, epoch).
///
/// Alias for [`OrbitState<Geocentric, GCRS>`][crate::astro::dynamics::OrbitState].
pub type CartesianState = crate::astro::dynamics::OrbitState<
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Julian Date epoch used by POD force models and propagators.
///
/// Alias for [`JulianDate`][crate::time::JulianDate].
pub type Epoch = crate::time::JulianDate;

// Convenience re-exports from the upstream canonical siderust types
/// Object-safe geocentric GCRS acceleration model trait object used by POD registries.
pub type DynSiderustForceModel = dyn crate::astro::dynamics::forces::AccelerationModel<
    crate::astro::dynamics::DynamicsContext,
    tempoch::TT,
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Concrete geocentric GCRS composite model type assembled from upstream accelerations.
pub type SiderustCompositeModel = crate::astro::dynamics::forces::CompositeModel<
    crate::astro::dynamics::DynamicsContext,
    tempoch::TT,
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Blanket trait for upstream accelerations defined on `DynamicsContext`, `TT`, geocentric center, and `GCRS`.
pub trait SiderustAccelerationModel:
    crate::astro::dynamics::forces::AccelerationModel<
        crate::astro::dynamics::DynamicsContext,
        tempoch::TT,
        crate::coordinates::centers::Geocentric,
        crate::coordinates::frames::GCRS,
    >
{
}

impl<T> SiderustAccelerationModel for T where
    T: crate::astro::dynamics::forces::AccelerationModel<
        crate::astro::dynamics::DynamicsContext,
        tempoch::TT,
        crate::coordinates::centers::Geocentric,
        crate::coordinates::frames::GCRS,
    >
{
}

#[doc(no_inline)]
pub use principia::{propagate_stm, StateTransitionMatrix};
#[doc(no_inline)]
pub use crate::astro::dynamics::{
    forces::{J2, TwoBody},
    DynamicsContext, OrbitState, Position, Propagator, PropagatorConfig, Velocity,
};
