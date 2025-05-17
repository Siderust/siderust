//! # Astro Module
//!
//! This module provides astronomical calculations and utilities, including:
//! - Sidereal time computations for tracking celestial objects.
//! - Nutation corrections for Earth's rotational axis wobble.
//! - Dynamical time calculations for precise timekeeping in astronomy.
//!
//! ## Submodules
//! - [`sidereal`]: Handles sidereal time calculations, including GST and LST.
//! - [`nutation`]: Computes nutation components (Δψ, Δε) and applies corrections.
//! - [`dynamical_time`]: Converts between Universal Time (UT) and Terrestrial Time (TT).
//! - [`orbit`]: Defines the `Orbit` struct for Keplerian orbital elements and computes heliocentric coordinates.

pub mod aberration;
pub mod proper_motion;
pub mod nutation;
pub mod precession;
pub mod sidereal;
pub mod dynamical_time;
pub mod orbit;
