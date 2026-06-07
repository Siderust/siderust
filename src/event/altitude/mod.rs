// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Unified Altitude Computation & Event API
//!
//! ## Scientific scope
//!
//! Topocentric altitude *h(t)* of a celestial body as observed from a
//! geodetic site on Earth, evaluated on the canonical TT axis. The body's
//! apparent direction is taken from the underlying analytical or numerical
//! engine (VSOP87/ELP2000/JPL ephemerides for solar-system bodies; ICRS
//! catalogue position for stars), so the regime of validity is bounded by
//! that of the chosen engine. Atmospheric refraction is **not** applied
//! here; thresholds such as `−0.833°` must be supplied by the caller.
//!
//! Event detection (rise/set crossings, culminations, time inside an
//! altitude band) uses Chebyshev-first crossing discovery with precise
//! validation and local scan+Brent fallback; temporal accuracy is limited
//! by the chosen `time_tolerance` and any explicit `scan_step_days`.
//!
//! ## Technical scope
//!
//! A clean, user‑friendly API for computing target altitude vs time and
//! finding events (crossings, culminations, altitude ranges) for **any**
//! celestial target.
//!
//! ## Module Structure
//!
//! - `types`, Core type definitions (events, queries, periods)
//! - `search`, Search options and configuration constants
//! - `compute`, Low-level altitude computation functions
//! - `events`, Event finding (crossings, culminations, ranges)
//! - `provider`, Trait-based dispatch for bodies
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
//! | `*_with_policy` variants | Same searches with an explicit apparent-position correction policy |
//!
//! All inputs/outputs are typed with `qtty` (`Degrees`, `JulianDate`, etc.).
//!
//! ## Time Scale
//!
//! `ModifiedJulianDate` / `Interval<ModifiedJulianDate>` values in this API are interpreted on
//! the TT axis. If your inputs are UTC timestamps, convert them with
//! `tempoch::Time::<tempoch::UTC>::from_chrono(...).to::<tempoch::TT>().into()`
//! into `ModifiedJulianDate` first.
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
//! use siderust::event::altitude::{crossings, AltitudePeriodsProvider, SearchOpts};
//! use siderust::bodies::Sun;
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::ECEF;
//! use siderust::time::{ModifiedJulianDate, Interval};
//! use siderust::qtty::*;
//!
//! let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
//! let window = Interval::new(
//!     siderust::ModifiedJulianDate::new(60000.0),
//!     siderust::ModifiedJulianDate::new(60001.0),
//! );
//!
//! // Pass any body that implements AltitudePeriodsProvider
//! let events = crossings(&Sun, &site, window, Degrees::new(0.0), SearchOpts::default());
//!
//! // Or use the trait methods directly
//! let alt_rad = Sun.altitude_at(&site, siderust::ModifiedJulianDate::new(60000.0));
//! ```
//!
//! ## References
//! None.

// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod events;
mod provider;
pub(crate) mod search;
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

pub use events::{
    above_threshold, above_threshold_with_policy, altitude_ranges, altitude_ranges_with_policy,
    below_threshold, below_threshold_with_policy, crossings, crossings_with_policy, culminations,
    culminations_with_policy,
};

// ---------------------------------------------------------------------------
// Re-exports: Trait & Provider Functions
// ---------------------------------------------------------------------------

pub use provider::{altitude_periods, AltitudePeriodsProvider};

#[doc(hidden)]
pub use events::{
    above_threshold_with_search_opts_v2, altitude_ranges_with_search_opts_v2,
    below_threshold_with_search_opts_v2, crossings_with_search_opts_v2,
};
