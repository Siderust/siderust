// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astro Module
//!
//! This module provides astronomical calculations and utilities, including:
//! - Sidereal time computations for tracking celestial objects.
//! - Nutation corrections for Earth's rotational axis wobble.
//! - Dynamical time calculations for precise timekeeping in astronomy.
//! - Annual Aberration to consider relativistic correction due to the Earth's orbital motion
//! - Precession and proper-motion helpers for mean/apparent coordinates.
//!
//! ## Submodules
//! - [`sidereal`]: Handles sidereal time calculations, including GST and LST.
//! - [`nutation`]: Computes nutation components (Δψ, Δε) and applies corrections.
//! - [`dynamical_time`]: Converts between Universal Time (UT) and Terrestrial Time (TT).
//! - [`aberration`]: Applies annual aberration corrections to celestial coordinates.
//! - [`orbit`]: Provides utilities for orbital mechanics and celestial object tracking.

pub mod aberration;
pub mod cio;
pub mod earth_rotation;
pub mod earth_rotation_provider;
pub mod eop;
pub mod era;
pub mod iers_data;
pub mod light_deflection;
pub mod nutation;
pub mod orbit;
pub mod polar_motion;
pub mod precession;
pub mod proper_motion;
pub mod sidereal;
