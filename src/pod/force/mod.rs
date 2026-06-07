// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # POD force-model configuration and registries
//!
//! ## Scientific scope
//! OD force-model configuration, registries, empirical accelerations, and thrust arcs.

pub mod config;
pub mod empirical_periodic;
pub mod low_thrust;
pub mod registry;
pub mod thrust;
pub mod thrust_arc;

pub use config::ForceModelConfig;
pub use empirical_periodic::{EmpiricalPeriodicAcceleration, PeriodicHarmonic};
pub use low_thrust::{LowThrustLog, LowThrustRecord};
pub use registry::{ForceModelFactory, ForceModelParams, ForceModelRegistry, ForceModelSpec};
pub use thrust::{mass_flow_rate, thrust_acceleration, ManeuverError, ThrustArc, G0_M_PER_S2};
pub use thrust_arc::ThrustArcConfig;

/// Geocentric GCRS Cartesian orbit state (position, velocity, epoch).
pub type CartesianState = crate::astro::dynamics::OrbitState<
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Julian Date epoch used by POD force models and propagators.
pub type Epoch = crate::time::JulianDate;

/// Object-safe geocentric GCRS acceleration model trait object.
pub type DynSiderustForceModel = dyn crate::astro::dynamics::forces::AccelerationModel<
    crate::astro::dynamics::DynamicsContext,
    tempoch::TT,
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Concrete geocentric GCRS composite model type.
pub type SiderustCompositeModel = crate::astro::dynamics::forces::CompositeModel<
    crate::astro::dynamics::DynamicsContext,
    tempoch::TT,
    crate::coordinates::centers::Geocentric,
    crate::coordinates::frames::GCRS,
>;

/// Blanket trait for upstream accelerations on DynamicsContext/TT/Geocentric/GCRS.
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
pub use crate::astro::dynamics::{
    forces::{TwoBody, J2},
    DynamicsContext, OrbitState, Position, Propagator, PropagatorConfig, Velocity,
};
