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
//! - Moon culmination detection (upper/lower meridian transits)
//!
//! ## Design Notes
//!
//! This module provides optimized lunar position routines that account for:
//! - **Topocentric parallax correction**: Critical for Moon due to its proximity (~1° at horizon)
//! - **Culmination-based search**: Partitions time into monotonic segments for efficient root-finding
//! - **Brent's method with endpoint reuse**: Avoids redundant ELP2000 evaluations
//!
//! ## Usage
//!
//! The main entry points are:
//! - [`moon_altitude_rad`]: Compute Moon altitude at a given Julian Date
//! - [`find_moon_altitude_periods_via_culminations`]: Find periods where Moon altitude satisfies a condition
//! - [`find_moon_above_horizon`] / [`find_moon_below_horizon`]: Convenience wrappers
//!
//! ---
//! _This module is under active development and will be expanded in future releases._

mod altitude_periods;

pub use altitude_periods::*;
