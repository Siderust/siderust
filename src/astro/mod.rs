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
//! - [`julian_date`]: Implements Julian Date calculations.
//! - [`modified_julian_date`]: Implements Modified Julian Date calculations.

pub mod aberration;
pub mod proper_motion;
pub mod nutation;
pub mod precession;
pub mod sidereal;
pub mod dynamical_time;
pub mod orbit;
pub mod julian_date;
pub mod modified_julian_date;

pub use julian_date::*;
pub use modified_julian_date::*;
