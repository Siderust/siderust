// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level propagation helpers built on `principia`.

pub mod propagator;

pub use principia::propagation::{
    EventDetector, EventOccurrence, PropagationError, PropagationResult, RadialThresholdEvent,
};
pub use propagator::{
    Dop853Propagator, Dopri5Propagator, FixedRk4Adapter, Propagator, PropagatorConfig,
    Rk4Propagator,
};
