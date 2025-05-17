//! # Time Units Module
//!
//! This module provides types and utilities for handling time-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! time systems and conversions between them.
//!
//! ## Features
//! - **Days**: A simple representation of time in days, with arithmetic operations.
//! - **Julian Day (JD)**: A continuous count of days since the beginning of the Julian Period (January 1, 4713 BCE).
//! - **Modified Julian Day (MJD)**: A variant of Julian Day used in technical and astronomical applications, defined as MJD = JD - 2400000.5.
//! - **Julian Year**: A standardized year of exactly 365.25 days, used in astronomy.
//! - **Centuries**: Representation of time intervals in Julian centuries (36525 days).
//! - **Years**: Representation of time intervals in years.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{JulianDay, ModifiedJulianDay, Days};
//!
//! // Create a Julian Day instance
//! let jd = JulianDay::new(2451545.0);
//!
//! // Convert to Modified Julian Day
//! let mjd = ModifiedJulianDay::new(jd.value() - 2400000.5);
//!
//! // Convert to UTC (if implemented)
//! let utc = jd.to_utc();
//!
//! // Perform arithmetic with Days
//! let future_jd = jd + Days::new(365.25); // Add one year
//! ```


mod days;
mod years;
mod julian_day;
mod julian_year;
mod modified_julian_day;
mod centuries;

pub use days::*;
pub use years::*;
pub use julian_day::*;
pub use julian_year::*;
pub use modified_julian_day::*;
pub use centuries::*;

