// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Result records produced by the propagation driver.

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
