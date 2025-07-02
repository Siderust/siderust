//! # Units Module
//!
//! This module provides a comprehensive set of strongly-typed units and utilities
//! for astronomical and scientific calculations. It is designed to ensure correctness,
//! clarity, and ease of use when working with various units of measurement.
//!
//! ## Features
//! - **Time Units**: Includes representations for Julian Day, Modified Julian Day, Days, Years, Julian Years, and Centuries.
//! - **Angular Units**: Provides types for Degrees, Radians, DMS (Degrees, Minutes, Seconds), HMS (Hours, Minutes, Seconds), and Arcseconds.
//! - **Length Units**: Includes types for meters and astronomical units (AU).
//! - **Velocity Units**: Provides types for meters per second and kilometers per second.
//! - **Mass Units**: Includes types for kilograms and solar masses.
//! - **Power Units**: Includes types for watts and solar luminosity.
//! - **Arithmetic Operations**: Supports arithmetic operations between compatible units, ensuring type safety.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{JulianDay, ModifiedJulianDay, Days, Degrees, Radians, DMS, Kilograms, SolarMass};
//!
//! // Time Units
//! let jd = JulianDay::new(2451545.0);
//! let mjd = ModifiedJulianDay::new(51544.5);
//! let days = Days::new(365.25);
//!
//! // Angular Units
//! let degrees = Degrees::new(180.0);
//! let radians = Radians::from(degrees);
//! let dms = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
//!
//! // Mass Units
//! let mass_kg = Kilograms::new(5.0);
//! let mass_solar = SolarMass::new(2.0);
//!
//! // Conversions
//! let jd_to_utc = jd.to_utc();
//! let mjd_to_jd = mjd.to_julian_day();
//! let dms_to_decimal = dms.as_f64();
//!
//! assert_eq!(radians.as_f64(), std::f64::consts::PI);
//! ```
//!
//! ## Modules
//! - [`time`]: Time-related units and utilities.
//! - [`angular`]: Angular measurement units and utilities.
//! - [`length`]: Length units and utilities.
//! - [`velocity`]: Velocity-related units and utilities.
//! - [`mass`]: Mass-related units and utilities.
//! - [`power`]: Power-related units and utilities.

pub(crate) mod arithmetic_ops;

pub mod angular;
pub mod time;
pub mod length;
pub mod velocity;
pub mod mass;
pub mod power;

pub use angular::*;
pub use time::*;
pub use length::*;
pub use velocity::*;
pub use mass::*;
pub use power::*;

use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign};
use simba::scalar::{ClosedAdd, ClosedSub};
use nalgebra::Scalar;

pub trait Unit:
    Copy
    + Clone
    + PartialEq
    + PartialOrd
    + std::fmt::Debug
    + std::fmt::Display
    + Add<Self>
    + Sub<Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + From<f64>
    + Into<f64>
    + PartialEq<f64>
    + AddAssign
    + SubAssign
    + ClosedAdd
    + ClosedSub
    + Scalar
    + From<AstronomicalUnit>
{
    const NAN: Self;

    fn sqrt(self) -> Self;
    fn powi(self, n: i32) -> Self;
    fn abs(self) -> Self;
}

impl Unit for f64 {
    const NAN: Self = f64::NAN;
    fn sqrt(self) -> Self { f64::sqrt(self) }
    fn powi(self, n: i32) -> Self { f64::powi(self, n) }
    fn abs(self) -> Self { f64::abs(self) }
}
