// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Unified Azimuth Computation & Event API
//!
//! A clean, user‑friendly API for computing target azimuth vs time and
//! finding events (bearing crossings, azimuth extrema, range periods) for
//! **any** celestial target.
//!
//! ## Module Structure
//!
//! - [`types`]    — Core type definitions (events, query)
//! - [`search`]   — Search options and configuration constants
//! - [`events`]   — Event finding (crossings, extrema, range periods)
//! - [`provider`] — Trait-based dispatch for bodies
//!
//! ## Public Functions
//!
//! | Function | Purpose |
//! |---|---|
//! | [`azimuth_crossings`] | When azimuth passes a specific bearing |
//! | [`azimuth_extrema`] | All local max/min (northernmost/southernmost) |
//! | [`azimuth_ranges`] | Intervals where azimuth is within `[min, max]` |
//! | [`in_azimuth_range`] | Same as `azimuth_ranges` (convenience alias) |
//! | [`outside_azimuth_range`] | Intervals where azimuth is outside `[min, max]` |
//!
//! All inputs/outputs are typed with `qtty` (`Degrees`, `ModifiedJulianDate`, etc.).
//!
//! ## Time Scale
//!
//! `ModifiedJulianDate` / `Period<MJD>` values in this API are interpreted on
//! the TT axis (`tempoch` canonical JD(TT) semantics).  If your inputs are UTC
//! timestamps, convert them with `ModifiedJulianDate::from_utc(…)` first.
//!
//! ## Azimuth Convention
//!
//! North-clockwise (North = 0°), values in `[0°, 360°)`.
//! A **wrap-around** query — e.g., 350° → 10° passing through North — is
//! expressed by setting `min_azimuth > max_azimuth` in [`AzimuthQuery`].
//!
//! ## Trait-Based API
//!
//! The [`AzimuthProvider`] trait provides a unified interface for computing
//! azimuth quantities for any body.  Implementations exist for
//! [`Sun`](crate::bodies::solar_system::Sun),
//! [`Moon`](crate::bodies::solar_system::Moon),
//! [`Star`](crate::bodies::Star), and
//! [`direction::ICRS`](crate::coordinates::spherical::direction::ICRS).
//!
//! ## Example
//!
//! ```rust
//! use siderust::calculus::azimuth::{azimuth_crossings, AzimuthProvider, SearchOpts};
//! use siderust::bodies::Sun;
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::ECEF;
//! use siderust::time::{ModifiedJulianDate, Period};
//! use qtty::*;
//!
//! let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
//! let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
//!
//! // Find when the Sun crosses due-South (180°):
//! let events = azimuth_crossings(&Sun, &site, window, Degrees::new(180.0), SearchOpts::default());
//! for e in &events {
//!     println!("Sun crosses South at MJD {:.6} ({:?})", e.mjd.value(), e.direction);
//! }
//!
//! // Find intervals where azimuth is between East (90°) and West (270°):
//! let query = siderust::calculus::azimuth::AzimuthQuery {
//!     observer: site,
//!     window,
//!     min_azimuth: Degrees::new(90.0),
//!     max_azimuth: Degrees::new(270.0),
//! };
//! let eastern_periods = Sun.azimuth_periods(&query);
//! ```

mod events;
mod provider;
mod search;
mod types;

// ---------------------------------------------------------------------------
// Re-exports: Types
// ---------------------------------------------------------------------------

pub use types::{
    AzimuthCrossingDirection, AzimuthCrossingEvent, AzimuthExtremum, AzimuthExtremumKind,
    AzimuthQuery,
};

// ---------------------------------------------------------------------------
// Re-exports: Search Options
// ---------------------------------------------------------------------------

pub use search::SearchOpts;

// ---------------------------------------------------------------------------
// Re-exports: Event Finding
// ---------------------------------------------------------------------------

pub use events::{
    azimuth_crossings, azimuth_extrema, azimuth_ranges, in_azimuth_range, outside_azimuth_range,
};

// ---------------------------------------------------------------------------
// Re-exports: Trait & Provider Functions
// ---------------------------------------------------------------------------

pub use provider::{azimuth_periods, AzimuthProvider};
