//! # Angular Units Module
//!
//! This module defines and handles various angular measurement units
//! used in astronomical calculations. It provides types and utilities
//! for converting between common angular formats such as:
//!
//! - **Degrees**: Decimal representation of an angle.
//! - **Radians**: Standard mathematical unit for angle measurement.
//! - **DMS (Degrees, Minutes, Seconds)**: Sexagesimal format used in many
//!   astronomical and geodetic contexts.
//! - **HMS (Hours, Minutes, Seconds)**: Commonly used in right ascension
//!   and time-based angular measurements.
//! - **Arcseconds**: Subdivision of degrees, used for high-precision measurements.
//!
//! ## Features
//! - Conversion between units (e.g., degrees to radians, DMS to degrees, HMS to degrees).
//! - Safe and precise parsing and formatting utilities.
//! - Strongly typed structs for each angular unit to avoid accidental misuse.
//! - Arithmetic operations for compatible angular units.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{Degrees, Radians, DMS, HMS};
//!
//! let deg = Degrees::new(180.0);
//! let rad = Radians::from(deg);
//! assert_eq!(rad.as_f64(), std::f64::consts::PI);
//!
//! let dms = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
//! let hms = HMS::new(5, 30, 0.0);
//! ```
//!
//! This module aims to make astronomical angle manipulations explicit,
//! correct, and ergonomic.


mod degrees;
mod radians;
mod dms;
mod hms;
mod arcsec;

pub use degrees::*;
pub use radians::*;
pub use dms::*;
pub use hms::*;
pub use arcsec::*;
