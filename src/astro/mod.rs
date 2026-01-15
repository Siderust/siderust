//! # Astro Module
//!
//! This module provides astronomical calculations and utilities, including:
//! - Sidereal time computations for tracking celestial objects.
//! - Nutation corrections for Earth's rotational axis wobble.
//! - Dynamical time calculations for precise timekeeping in astronomy.
//! - Annual Aberration to consider relativistic correction due to the Earth's orbital motion
//! - Astronomical Dates such as JulianDate and Modified Julian Date.
//!
//! ## Submodules
//! - [`sidereal`]: Handles sidereal time calculations, including GST and LST.
//! - [`nutation`]: Computes nutation components (Δψ, Δε) and applies corrections.
//! - [`dynamical_time`]: Converts between Universal Time (UT) and Terrestrial Time (TT).
//! - [`aberration`]: Applies annual aberration corrections to celestial coordinates.
//! - [`orbit`]: Provides utilities for orbital mechanics and celestial object tracking.

pub mod aberration;
pub mod dynamical_time;
pub mod nutation;
pub mod orbit;
pub mod precession;
pub mod proper_motion;
pub mod sidereal;

// Re-export time types for backward compatibility
pub use crate::time::{JulianDate, ModifiedJulianDate};

