// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Calculus Module
//!
//! This module provides algorithms and data structures for calculations related to the Moon,
//! focusing on the Moon's apparent position as seen from Earth.
//!
//! ## Current Features
//!
//! - Calculation of the Moon's geocentric and topocentric coordinates using ELP2000-82B theory
//! - Optimized moon altitude calculations with topocentric parallax correction
//! - Finding periods when Moon is above/below specific altitude thresholds
//! - Altitude range detection (`find_moon_altitude_range`)
//!
//! ## Design Notes
//!
//! All period-finding delegates to [`crate::calculus::math_core::intervals`]
//! which provides scan + Brent refinement + interval assembly.  This module
//! supplies the Moon-altitude closure and JD↔Mjd conversions.
//!
//! ## Usage
//!
//! The main entry points are:
//! - [`moon_altitude_rad`]: Compute Moon altitude at a given Julian Date
//! - [`find_moon_above_horizon`] / [`find_moon_below_horizon`]: Above/below threshold
//! - [`find_moon_altitude_range`]: Within a min/max altitude band
//!
//! ---
//! _This module is under active development and will be expanded in future releases._

mod altitude_periods;
pub(crate) mod moon_cache;
mod moon_equations;

pub use altitude_periods::*;
pub use moon_cache::MoonAltitudeContext;
