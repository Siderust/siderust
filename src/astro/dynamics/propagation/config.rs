// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Low-level driver configuration for the [`super::driver::propagate`] function.
//!
//! ## Scope
//!
//! Provides [`PropagationConfig<C, F>`], the typed configuration bundle for the
//! adaptive propagation driver: time bounds, step hints, tolerances, event
//! detectors, output cadence, and step budgeting.
//!
//! See [`super::propagator::Propagator`] and the high-level API for the
//! user-facing wrapper.
//!
//! ## Typical usage
//!
//! ```rust,ignore
//! use siderust::astro::dynamics::propagation::PropagationConfig;
//! use siderust::time::{Time, TT, JulianDate};
//! use siderust::qtty::Second;
//!
//! let t_start = Time::<TT>::from_jd(JulianDate::new(2_451_545.0));
//! let t_end = t_start + Second::new(86_400.0);
//!
//! let cfg = PropagationConfig::new(t_start, t_end)
//!     .with_initial_step(Second::new(30.0))
//!     .with_max_step(Second::new(600.0))
//!     .with_output_every(Second::new(60.0));
//! ```

use super::events::EventDetector;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::ext_qtty::tolerances::IntegratorTolerances;
use crate::qtty::Second;
use crate::time::{Time, TT};

/// Low-level driver configuration.  Use [`super::propagator::PropagatorConfig`]
/// + [`super::propagator::Propagator`] for the ergonomic high-level API.
pub struct PropagationConfig<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub t_start: Time<TT>,
    pub t_end: Time<TT>,
    pub h0: Second,
    pub h_max: Second,
    pub h_min: Second,
    pub tolerances: IntegratorTolerances,
    pub output_every: Option<Second>,
    pub output_at: Vec<Time<TT>>,
    pub events: Vec<Box<dyn EventDetector<C, F>>>,
    pub max_steps: u32,
}

impl<C, F> PropagationConfig<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn new(t_start: Time<TT>, t_end: Time<TT>) -> Self {
        Self {
            t_start,
            t_end,
            h0: Second::new(30.0),
            h_max: Second::new(86_400.0),
            h_min: Second::new(1e-6),
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-3, 1e-6),
            output_every: None,
            output_at: Vec::new(),
            events: Vec::new(),
            max_steps: 1_000_000,
        }
    }

    pub fn with_initial_step(mut self, h0: Second) -> Self {
        self.h0 = h0;
        self
    }
    pub fn with_max_step(mut self, h_max: Second) -> Self {
        self.h_max = h_max;
        self
    }
    pub fn with_min_step(mut self, h_min: Second) -> Self {
        self.h_min = h_min;
        self
    }
    pub fn with_tolerances(mut self, tols: IntegratorTolerances) -> Self {
        self.tolerances = tols;
        self
    }
    pub fn with_output_every(mut self, dt: Second) -> Self {
        self.output_every = Some(dt);
        self
    }
    pub fn with_output_at(mut self, times: Vec<Time<TT>>) -> Self {
        self.output_at = times;
        self
    }
    pub fn with_event(mut self, ev: Box<dyn EventDetector<C, F>>) -> Self {
        self.events.push(ev);
        self
    }
    pub fn with_max_steps(mut self, max_steps: u32) -> Self {
        self.max_steps = max_steps;
        self
    }

    /// Returns the total signed duration `t_end - t_start` in seconds.
    pub fn total_duration_s(&self) -> f64 {
        (self.t_end - self.t_start).value()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::propagation::events::AltitudeEvent;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::ext_qtty::tolerances::IntegratorTolerances;
    use crate::qtty::{Kilometers, Second};
    use crate::time::{JulianDate, Time, TT};

    fn epoch(jd: f64) -> Time<TT> {
        let s = OrbitState::new_at_jd(
            JulianDate::new(jd),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        s.epoch
    }

    #[test]
    fn new_sets_defaults() {
        let cfg: PropagationConfig = PropagationConfig::new(epoch(2_451_545.0), epoch(2_451_546.0));
        assert!((cfg.h0.value() - 30.0).abs() < 1e-10);
        assert!(cfg.output_every.is_none());
        assert!(cfg.events.is_empty());
        assert_eq!(cfg.max_steps, 1_000_000);
    }

    #[test]
    fn builder_methods_set_fields() {
        let cfg: PropagationConfig = PropagationConfig::new(epoch(2_451_545.0), epoch(2_451_546.0))
            .with_initial_step(Second::new(10.0))
            .with_max_step(Second::new(300.0))
            .with_min_step(Second::new(0.001))
            .with_tolerances(IntegratorTolerances::uniform(1e-8, 1e-5, 1e-8))
            .with_output_every(Second::new(60.0))
            .with_max_steps(50_000);
        assert!((cfg.h0.value() - 10.0).abs() < 1e-10);
        assert!((cfg.h_max.value() - 300.0).abs() < 1e-10);
        assert!((cfg.h_min.value() - 0.001).abs() < 1e-14);
        assert!(cfg.output_every.is_some());
        assert!((cfg.output_every.unwrap().value() - 60.0).abs() < 1e-10);
        assert_eq!(cfg.max_steps, 50_000);
    }

    #[test]
    fn with_output_at_sets_vec() {
        let times = vec![epoch(2_451_545.01), epoch(2_451_545.02)];
        let cfg: PropagationConfig =
            PropagationConfig::new(epoch(2_451_545.0), epoch(2_451_545.1)).with_output_at(times);
        assert_eq!(cfg.output_at.len(), 2);
    }

    #[test]
    fn with_event_appends_event() {
        let event = AltitudeEvent::new(Kilometers::new(600.0), Kilometers::new(6_371.0));
        let cfg: PropagationConfig = PropagationConfig::new(epoch(2_451_545.0), epoch(2_451_546.0))
            .with_event(Box::new(event));
        assert_eq!(cfg.events.len(), 1);
    }

    #[test]
    fn total_duration_s_positive() {
        let dt_days = 1.0;
        let dt_s = dt_days * 86_400.0;
        let cfg: PropagationConfig =
            PropagationConfig::new(epoch(2_451_545.0), epoch(2_451_545.0 + dt_days));
        let dur = cfg.total_duration_s();
        assert!((dur - dt_s).abs() < 1.0, "total duration off: {dur}");
    }

    #[test]
    fn total_duration_s_negative_for_backward() {
        let cfg: PropagationConfig = PropagationConfig::new(epoch(2_451_546.0), epoch(2_451_545.0));
        assert!(
            cfg.total_duration_s() < 0.0,
            "backward propagation must give negative duration"
        );
    }
}
