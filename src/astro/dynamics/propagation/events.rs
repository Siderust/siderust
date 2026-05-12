// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Zero-crossing event detectors for the propagation driver.
//!
//! ## Scope
//!
//! Provides the [`EventDetector<C, F>`] trait and built-in implementations
//! for altitude-based event detection (e.g., apogee, perigee, eclipse).
//!
//! ## Equations
//!
//! Event detection uses a zero-crossing function `g(t, x)`:
//! - When `g(t) changes sign between consecutive steps, an event occurred.
//! - If `terminal()` returns `true`, propagation terminates immediately.
//!
//! Example (altitude trigger):
//! ```text
//! g(t) = |r(t)| − (R_body + h_trigger)
//! ```
//! - `g > 0` means above trigger altitude.
//! - `g < 0` means below trigger altitude.
//!
//! ## Built-in detectors
//!
//! | Type | Switching function |
//! |------|-------------------|
//! | [`AltitudeEvent`] | `\|r\| − (R + h)` |
//!
//! ## Units & frames
//!
//! All positions in km (geocentric).  Altitudes in km above a spherical body.
//!
//! ## Custom detectors
//!
//! Implement [`EventDetector`] to add custom zero-crossing functions.

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::Kilometers;

/// Zero-crossing detector for a scalar switching function.
pub trait EventDetector<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn name(&self) -> &'static str;
    /// Switching function `g(t, x)`.  Detection occurs at sign changes.
    fn evaluate(
        &self,
        state: &OrbitState<C, F>,
        ctx: &DynamicsContext,
    ) -> Result<f64, DynamicsError>;
    /// If `true`, propagation terminates immediately when this event fires.
    fn terminal(&self) -> bool {
        false
    }
}

/// Detects when the spacecraft crosses a specified altitude above a spherical body.
///
/// The switching function is `|r| - (body_radius + trigger_altitude)`.
/// A positive value means the spacecraft is above the trigger altitude.
pub struct AltitudeEvent {
    pub trigger_altitude: Kilometers,
    pub body_radius: Kilometers,
    pub terminal: bool,
}

impl AltitudeEvent {
    pub fn new(trigger_altitude: Kilometers, body_radius: Kilometers) -> Self {
        Self {
            trigger_altitude,
            body_radius,
            terminal: false,
        }
    }
    pub fn terminal(mut self, terminal: bool) -> Self {
        self.terminal = terminal;
        self
    }
}

impl EventDetector<Geocentric, GCRS> for AltitudeEvent {
    fn name(&self) -> &'static str {
        "altitude"
    }
    fn evaluate(
        &self,
        state: &OrbitState<Geocentric, GCRS>,
        _ctx: &DynamicsContext,
    ) -> Result<f64, DynamicsError> {
        let x = state.position.x().value();
        let y = state.position.y().value();
        let z = state.position.z().value();
        let r = (x * x + y * y + z * z).sqrt();
        Ok(r - (self.body_radius.value() + self.trigger_altitude.value()))
    }
    fn terminal(&self) -> bool {
        self.terminal
    }
}
