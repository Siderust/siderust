// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Calculus Module
//!
//! ## Scientific scope
//!
//! Algorithms and data structures for the apparent position of the Sun as
//! seen from a topocentric Earth observer. Covers altitude (no
//! refraction), azimuth, twilight phase classification (civil, nautical,
//! astronomical), and time intervals during which the Sun is above/below a
//! threshold or inside a band. Validity is bounded by the underlying
//! VSOP87 model used for the Sun–Earth geometry.
//!
//! ## Technical scope
//!
//! - Apparent geocentric equatorial coordinates with nutation + aberration.
//! - Optimised Sun altitude/azimuth closures and internal period helpers
//!   (`solar_*_impl`, `sun_altitude_rad`) built on labelled crossings.
//! - Twilight phase classification via [`twilight_classification`].
//!
//! All period‑finding delegates to `crate::event::search::intervals`
//! which provides scan + Brent refinement + interval assembly. This module
//! supplies the Sun‑altitude closure and JD↔MJD conversions.
//!
//! ## References
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann‑Bell.
//! - Bretagnon, P. & Francou, G. (1988). "Planetary theories in
//!   rectangular and spherical variables. VSOP 87 solutions". *A&A*,
//!   202, 309–315.

mod daily_events;
mod sun_equations;

pub(crate) mod altitude_periods;
pub(crate) mod azimuth;
pub mod classification;
pub mod night_types;

pub(crate) use altitude_periods::*;
pub(crate) use azimuth::*;
pub use classification::{twilight_classification, TwilightPhase};
pub use night_types::*;
