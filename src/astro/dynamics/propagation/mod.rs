// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level adaptive propagation driver with event detection.
//!
//! ## Scope
//!
//! Provides the [`propagate`] function, which composes any [`AdaptiveStepper`]
//! with a [`PropagatorConfig`] to handle event detection, output sampling,
//! and step budgeting in a single high-level call.
//!
//! ## Workflow
//!
//! 1. Create a [`PropagatorConfig`] with time bounds, tolerances, event detectors.
//! 2. Call [`propagate`] with an integrator (DOPRI5 or DOP853), force model, and context.
//! 3. Inspect the [`PropagationResult`] for samples, event times, and step counts.
//!
//! ## Failure modes
//!
//! - [`PropagationError::StepBelowMinimum`] if the integrator shrinks below `h_min`.
//! - [`PropagationError::MaxStepsExceeded`] if the step budget is exhausted.
//! - [`PropagationError::EventEvaluation`] if an event detector's switching function fails.
//! - [`PropagationError::StepControl`] for any integrator error.
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §4.

pub mod config;
pub mod driver;
pub mod error;
pub mod events;
pub mod propagator;
pub mod result;

pub use config::PropagationConfig;
pub use driver::propagate;
pub use error::PropagationError;
pub use events::{AltitudeEvent, EventDetector};
pub use propagator::{IntegratorChoice, Propagator, PropagatorConfig};
pub use result::{EventOccurrence, PropagationResult};
