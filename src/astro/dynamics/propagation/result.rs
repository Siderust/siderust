// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Result records produced by the propagation driver.
//!
//! ## Scope
//!
//! Provides [`PropagationResult<C, F>`] — the output structure from
//! [`super::driver::propagate`], containing propagated samples, event
//! occurrences, and step statistics.
//!
//! ## Structure
//!
//! ```text
//! PropagationResult {
//!     samples: Vec<OrbitState<C, F>>,       // Initial + output states
//!     events: Vec<EventOccurrence<C, F>>,   // Event times and states
//!     steps_taken: u32,                     // Number of accepted steps
//!     steps_rejected: u32,                  // Number of rejected steps
//! }
//! ```
//!
//! The `samples` vector always includes the initial state (at `t_start`).
//! Additional entries are added if `output_every` or `output_at` was configured.
//!
//! ## Usage
//!
//! ```rust,ignore
//! let result = propagate(&integrator, &force, initial_state, &cfg, &ctx)?;
//! for sample in &result.samples {
//!     println!("{:?}", sample.position);
//! }
//! for event in &result.events {
//!     println!("Event '{}' at {:?}", event.event_name, event.state.epoch);
//! }
//! ```

use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};

#[derive(Debug, Clone)]
pub struct EventOccurrence<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub event_name: &'static str,
    pub state: OrbitState<C, F>,
}

#[derive(Debug, Clone)]
pub struct PropagationResult<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Final state and any output samples requested via the [`super::config::PropagatorConfig`].
    pub samples: Vec<OrbitState<C, F>>,
    pub events: Vec<EventOccurrence<C, F>>,
    pub steps_taken: u32,
    pub steps_rejected: u32,
}
