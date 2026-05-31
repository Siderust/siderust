// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Runtime mission context
//!
//! ## Scientific scope
//!
//! The [`MissionContext`] type aggregates the typed objects an
//! operational pipeline needs at runtime: the catalog of ground-station
//! [`crate::mission::site::Location`]s and the catalog of on-board
//! [`crate::mission::geometry::Instrument`]s. It is intentionally a passive container of
//! typed handles — there is no service registry, no async runtime, no
//! network layer. Those operational concerns live in **SatOps**.
//!
//! ## Technical scope
//!
//! - Build a context with [`MissionContext::new`] (no items registered)
//!   or [`MissionContext::default`] and progressively call the `add_*` /
//!   `with_*` builder methods.
//! - Look up locations / instruments by string alias; aliases are
//!   case-sensitive and unique per category.
//! - The context owns its handles; cloning is cheap because the
//!   contained types are themselves cheaply clonable.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §1.2 (mission analysis context).

use std::collections::HashMap;

use crate::mission::geometry::Instrument;
use crate::mission::site::Location;

/// Runtime mission-analysis context.
#[derive(Debug, Default, Clone)]
pub struct MissionContext {
    locations: HashMap<String, Location>,
    instruments: HashMap<String, Instrument>,
}

impl MissionContext {
    /// Create an empty mission context.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a location under a free-text alias.
    pub fn add_location(&mut self, alias: impl Into<String>, location: Location) {
        self.locations.insert(alias.into(), location);
    }

    /// Register an instrument under a free-text alias.
    pub fn add_instrument(&mut self, alias: impl Into<String>, instrument: Instrument) {
        self.instruments.insert(alias.into(), instrument);
    }

    /// Builder variant of [`Self::add_location`].
    pub fn with_location(mut self, alias: impl Into<String>, location: Location) -> Self {
        self.add_location(alias, location);
        self
    }

    /// Builder variant of [`Self::add_instrument`].
    pub fn with_instrument(mut self, alias: impl Into<String>, instrument: Instrument) -> Self {
        self.add_instrument(alias, instrument);
        self
    }

    /// Look up a registered location by alias.
    pub fn location(&self, alias: &str) -> Option<&Location> {
        self.locations.get(alias)
    }

    /// Look up a registered instrument by alias.
    pub fn instrument(&self, alias: &str) -> Option<&Instrument> {
        self.instruments.get(alias)
    }

    /// Remove the location registered under `alias`, returning it
    /// if present.
    pub fn remove_location(&mut self, alias: &str) -> Option<Location> {
        self.locations.remove(alias)
    }

    /// Remove the instrument registered under `alias`, returning it
    /// if present.
    pub fn remove_instrument(&mut self, alias: &str) -> Option<Instrument> {
        self.instruments.remove(alias)
    }

    /// Iterate over `(alias, location)` pairs.
    pub fn locations(&self) -> impl Iterator<Item = (&String, &Location)> {
        self.locations.iter()
    }

    /// Iterate over `(alias, instrument)` pairs.
    pub fn instruments(&self) -> impl Iterator<Item = (&String, &Instrument)> {
        self.instruments.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mission::geometry::{Fov, LocalFrame};
    use qtty::Quantity;

    #[test]
    fn empty_context_has_no_aliases() {
        let ctx = MissionContext::new();
        assert!(ctx.location("foo").is_none());
        assert!(ctx.instrument("bar").is_none());
    }

    #[test]
    fn add_and_lookup_round_trip() {
        let frame = LocalFrame::new(Quantity::new(0.5), Quantity::new(0.3));
        let loc = Location::new("VLBI-A", frame);
        let inst = Instrument::new(
            "cam",
            Fov::Conical {
                half_angle: Quantity::new(0.1),
            },
        );
        let ctx = MissionContext::new()
            .with_location("station-a", loc.clone())
            .with_instrument("cam-1", inst.clone());
        assert_eq!(ctx.location("station-a").unwrap().id, "VLBI-A");
        assert_eq!(ctx.instrument("cam-1").unwrap().id, "cam");
    }

    #[test]
    fn remove_returns_registered_object() {
        let frame = LocalFrame::new(Quantity::new(0.5), Quantity::new(0.3));
        let loc = Location::new("VLBI-A", frame);
        let mut ctx = MissionContext::new();
        ctx.add_location("a", loc);
        assert!(ctx.remove_location("a").is_some());
        assert!(ctx.remove_location("a").is_none());
    }
}
