// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Spacecraft dynamics
//!
//! Cartesian state propagation, force models, numerical integrators, and
//! ancillary science models for satellite orbit mechanics.
//!
//! ## Scientific scope
//!
//! While [`crate::astro::orbit`] models orbits as Keplerian element sets,
//! many high-fidelity applications — orbit propagation under non-Keplerian
//! perturbations (J2, third-body, drag, SRP), state-transition matrices,
//! covariance transport, precise orbit determination — work directly with a
//! Cartesian inertial state `(r, v)`. This module provides the typed primitive
//! for that representation, [`OrbitState`], plus the surrounding spacecraft
//! description (mass, drag and SRP cross-sections).
//!
//! ## Technical scope
//!
//! All position, velocity, and acceleration values are [`affn`] vectors with
//! [`qtty`] units, so frame and unit constraints are checked at compile time.
//! Time derivatives are carried as the typed [`StateDerivative`] aggregate.
//!
//! Submodules host:
//!
//! - [`forces`] — [`ForceModel`] trait and standard models (two-body, J2, drag, SRP, third-body, geopotential, relativity).
//! - [`integrators`] — RK4 (fixed-step) and DOPRI5 / DOP853 (adaptive).
//! - [`propagation`] — high-level adaptive driver with event detection.
//! - [`variational`] — analytic state-transition-matrix propagation via variational equations.
//! - [`gravity`] — geopotential field provider abstraction and acceleration kernel.
//! - [`atmosphere`] — atmospheric density provider abstraction.
//! - [`frames`] — local orbital frames (RTN, VNC, LVLH).
//! - [`covariance`] — frame-tagged 6×6 state covariance with transport.
//! - [`stm`] — finite-difference state-transition matrix computation.
//!
//! ## Units and frames
//!
//! | Quantity | Type/Units | Convention |
//! |----------|-----------|-----------|
//! | Position | km | Geocentric, GCRS |
//! | Velocity | km/s | GCRS |
//! | Acceleration | km/s² | GCRS (unless stated otherwise) |
//! | GM | km³/s² | EGM2008/WGS-84 |
//! | Density | kg/m³ | SI (geodetic altitude from sea level) |
//! | Time | seconds (TT) | Typed as [`qtty::Second`] |
//!
//! ## Failure modes
//!
//! Most errors are surface through [`DynamicsError`]:
//! - Provider unavailability (ephemeris, EOP, gravity, atmosphere)
//! - Degenerate geometry (zero position, zero velocity, r ∥ v)
//! - Geopotential degree/order out of range
//! - Integration step control failures
//!
//! See [`errors`] for details.

pub mod atmosphere;
pub mod context;
pub mod covariance;
pub mod errors;
pub mod forces;
pub mod frames;
pub mod gravity;
pub mod integrators;
pub mod propagation;
pub mod state;
pub mod stm;
pub mod units;
pub mod variational;

pub use context::{
    Conventions, DynamicsContext, DynamicsContextBuilder, EarthOrientationProvider,
    PrecessionModel, SolarActivityProvider,
};
pub use covariance::ProcessNoise;
pub use errors::{DynamicsError, LocalFrameError};
pub use state::{
    Acceleration, AccelerationUnit, OrbitState, Position, SpacecraftProperties, SpacecraftState,
    StateDerivative, Velocity, VelocityUnit,
};
#[allow(deprecated)]
pub use stm::{finite_diff_stm, finite_diff_stm_series, StateTransition};
pub use units::GravitationalParameter;

pub use propagation::{
    AltitudeEvent, EventDetector, EventOccurrence, PropagationConfig, PropagationResult,
    Propagator, PropagatorConfig, IntegratorChoice,
};
pub use variational::{propagate_stm, StateTransitionMatrix};
