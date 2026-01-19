//! # Solar Calculus Module
//!
//! This module provides algorithms and data structures for calculations related to the Sun
//! and solar system barycenter, focusing on the Sun's apparent position as seen from Earth.
//!
//! ## Current Features
//!
//! - Calculation of the Sun's apparent geocentric equatorial coordinates, including
//!   corrections for nutation and aberration.
//!
//! ## Design Notes
//!
//! This module is intentionally minimal in its current form, providing only the most
//! essential solar position routines. It is designed for extensibility and will be expanded
//! in future releases to include:
//! - Solar elongation and phase angle calculations
//! - Solar limb darkening and apparent diameter
//! - Solar ephemerides for arbitrary epochs
//! - Solar system barycenter computations
//! - Additional solar-related phenomena and corrections
//!
//! ## Usage
//!
//! The main entry point is the [`crate::bodies::solar_system::Sun`] type, which exposes methods for obtaining the Sun's
//! apparent position `get_apparent_geocentric_equ`.
//!
//! ---
//! _This module is under active development and will be expanded in future releases._

mod sun_equations;

pub mod altitude_periods;

pub use altitude_periods::*;
