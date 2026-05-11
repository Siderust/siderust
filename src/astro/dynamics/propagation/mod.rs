// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level adaptive propagation driver and supporting types.

pub mod config;
pub mod driver;
pub mod error;
pub mod events;
pub mod result;

pub use config::PropagatorConfig;
pub use driver::propagate;
pub use error::PropagationError;
pub use events::{AltitudeEvent, EventDetector};
pub use result::{EventOccurrence, PropagationResult};
