// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Calculus Module
//!
//! This module provides algorithms and data structures for calculations related to the Sun
//! and solar system barycenter, focusing on the Sun's apparent position as seen from Earth.
//!
//! ## Current Features
//!
//! - Calculation of the Sun's apparent geocentric equatorial coordinates, including
//!   corrections for nutation and aberration.
//! - Optimized sun altitude calculations
//! - Finding periods when Sun is above/below specific altitude thresholds (day/night)
//! - Altitude range detection (`find_sun_range_periods`)
//!
//! ## Design Notes
//!
//! All period-finding delegates to [`crate::calculus::math_core::intervals`]
//! which provides scan + Brent refinement + interval assembly.  This module
//! supplies the Sun-altitude closure and JD↔Mjd conversions.
//!
//! ## Usage
//!
//! The main entry points are:
//! - [`sun_altitude_rad`]: Compute Sun altitude at a given Julian Date
//! - [`find_day_periods`] / [`find_night_periods`]: Above/below threshold
//! - [`find_sun_range_periods`]: Within a min/max altitude band
//!
//! ---
//! _This module is under active development and will be expanded in future releases._

mod sun_equations;

pub(crate) mod altitude_periods;
pub mod night_types;

pub(crate) use altitude_periods::*;
pub use night_types::*;
