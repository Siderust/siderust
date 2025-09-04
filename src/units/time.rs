//! # Time Units Module
//!
//! This module provides types and utilities for handling time-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! time systems and conversions between them.
//!
//! ## Features
//! - **Days**: A simple representation of time in days, with arithmetic operations.
//! - **Julian Year**: A standardized year of exactly 365.25 days, used in astronomy.
//! - **Centuries**: Representation of time intervals in Julian centuries (36525 days).
//! - **Years**: Representation of time intervals in years.

use super::*;

pub enum Time {}
impl Dimension for Time {}

pub trait TimeUnit: Unit<Dim = Time> {}
impl<T: Unit<Dim = Time>> TimeUnit for T {}

define_unit!("ms", Millisecond, Time, 1.0 / (24.0 * 3600.0 * 1_000.0));
pub type Milliseconds = Quantity<Millisecond>;
pub const MILLISEC: Milliseconds = Milliseconds::new(1.0);

define_unit!("sec", Second, Time, 1.0 / (24.0 * 3600.0));
pub type Seconds = Quantity<Second>;
pub const SEC: Seconds = Seconds::new(1.0);

define_unit!("min", Minute, Time, 1.0 / (24.0 * 60.0));
pub type Minutes = Quantity<Minute>;
pub const MIN: Minutes = Minutes::new(1.0);

define_unit!("h", Hour, Time, 1.0 / 24.0);
pub type Hours = Quantity<Hour>;
pub const HOUR: Hours = Hours::new(1.0);

// Mean solar day
define_unit!("d", Day, Time, 1.0);
pub type Days = Quantity<Day>;
pub const DAY: Days = Days::new(1.0);

// Week: 7 solar
define_unit!("wk", Week, Time, 7.0);
pub type Weeks = Quantity<Week>;
pub const WEEK: Weeks = Weeks::new(1.0);

// Mean tropical year (IAU 2015)
define_unit!("yr", Year, Time, 365.242_5);
pub type Years = Quantity<Year>;
pub const YEAR: Years = Years::new(1.0);

// Century: 100 mean tropical years
define_unit!("cent", Century, Time, 36_524.25);
pub type Centuries = Quantity<Century>;
pub const CENTURY: Centuries = Centuries::new(1.0);

define_unit!("JY", JulianYear, Time, 365.25);
pub type JulianYears = Quantity<JulianYear>;
pub const JULIAN_YEAR: JulianYears = JulianYears::new(1.0);

// Julian century: exactly 36,525 days (365.25 × 100).
define_unit!("JC", JulianCentury, Time, 36_525.0);
pub type JulianCenturies = Quantity<JulianCentury>;
pub const JULIAN_CENTURY: JulianCenturies = JulianCenturies::new(1.0);
