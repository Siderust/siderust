// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Unified Altitude Computation & Event API
//!
//! A clean, user‑friendly API for computing target altitude vs time and
//! finding events (crossings, culminations, altitude ranges) for **any**
//! celestial target.
//!
//! ## Module Structure
//!
//! - [`types`] — Core type definitions (events, queries, periods)
//! - [`search`] — Search options and configuration constants
//! - [`compute`] — Low-level altitude computation functions
//! - [`events`] — Event finding (crossings, culminations, ranges)
//! - [`provider`] — Trait-based dispatch for bodies
//!
//! ## Public Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`crossings`] | Threshold crossings in a time window |
//! | [`culminations`] | All local maxima/minima of altitude |
//! | [`altitude_ranges`] | Intervals where altitude is within `[h_min, h_max]` |
//! | [`above_threshold`] | Intervals where altitude is above threshold |
//! | [`below_threshold`] | Intervals where altitude is below threshold |
//!
//! All inputs/outputs are typed with `qtty` (`Degrees`, `JulianDate`, etc.).
//!
//! ## Time Scale
//!
//! All `ModifiedJulianDate` / `Period<MJD>` values in this module are on
//! tempoch's canonical JD(TT) axis (not raw UTC-MJD counters).
//! Convert UTC timestamps with `ModifiedJulianDate::from_utc(...)`
//! before calling these APIs.
//!
//! ## Trait-Based API
//!
//! The [`AltitudePeriodsProvider`] trait provides a unified interface for
//! computing altitude periods of any celestial body. Implementations exist
//! for [`Sun`](crate::bodies::solar_system::Sun),
//! [`Moon`](crate::bodies::solar_system::Moon),
//! [`Star`](crate::bodies::Star), and
//! [`direction::ICRS`](crate::coordinates::spherical::direction::ICRS).
//!
//! All event-finding functions are generic over `AltitudePeriodsProvider`,
//! so you can pass any supported body directly.
//!
//! ## Example
//!
//! ```rust
//! use siderust::calculus::altitude::{crossings, AltitudePeriodsProvider, SearchOpts};
//! use siderust::bodies::Sun;
//! use siderust::coordinates::centers::ObserverSite;
//! use siderust::time::{ModifiedJulianDate, Period};
//! use qtty::*;
//!
//! let site = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
//! let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
//!
//! // Pass any body that implements AltitudePeriodsProvider
//! let events = crossings(&Sun, &site, window, Degrees::new(0.0), SearchOpts::default());
//!
//! // Or use the trait methods directly
//! let alt_rad = Sun.altitude_at(&site, siderust::time::ModifiedJulianDate::new(60000.0));
//! ```

// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod events;
mod provider;
mod search;
mod types;

// ---------------------------------------------------------------------------
// Re-exports: Types
// ---------------------------------------------------------------------------

pub use types::{
    AltitudeQuery, CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind,
};

// ---------------------------------------------------------------------------
// Re-exports: Search Options
// ---------------------------------------------------------------------------

pub use search::SearchOpts;

// ---------------------------------------------------------------------------
// Re-exports: Core Computation
// ---------------------------------------------------------------------------

// fixed_star_altitude_rad moved to stellar module

// ---------------------------------------------------------------------------
// Re-exports: Event Finding
// ---------------------------------------------------------------------------

pub use events::{above_threshold, altitude_ranges, below_threshold, crossings, culminations};

// ---------------------------------------------------------------------------
// Re-exports: Trait & Provider Functions
// ---------------------------------------------------------------------------

pub use provider::{altitude_periods, AltitudePeriodsProvider};
