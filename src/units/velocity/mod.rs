//! # Velocity Units Module
//!
//! This module provides types and utilities for handling velocity-related calculations
//! in astronomical and scientific contexts. It includes representations for angular
//! velocity such as radians per day and conversions from angular and time units.
//!
//! ## Features
//! - **RadiansPerDay**: Represents angular velocity in radians per day.
//! - Arithmetic operations for velocity types.
//! - Conversion from `Radians` and `Days`.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{Radians, Days, RadiansPerDay};
//!
//! let angle = Radians::new(std::f64::consts::PI);
//! let time = Days::new(2.0);
//! let velocity = angle / time; // RadiansPerDay
//! assert_eq!(velocity.value(), std::f64::consts::PI / 2.0);
//! ```

mod radians_per_day;
mod dms_per_year;
mod au_per_day;

pub use radians_per_day::*;
pub use dms_per_year::*;
pub use au_per_day::*;
