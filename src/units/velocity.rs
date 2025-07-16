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

use crate::units::*;

pub enum Velocity {}
impl Dimension for Velocity {}
pub trait VelocityUnit: Unit<Dim = Velocity> {}
impl<T: Unit<Dim = Velocity>> VelocityUnit for T {}


define_unit!("m/s", MeterPerSec, Velocity, Meter::RATIO/Second::RATIO);
pub type MetersPerSec = Quantity<MeterPerSec>;

define_unit!("Km/s", KilometerPerSec, Velocity, Kilometer::RATIO/Second::RATIO);
pub type KilometersPerSec = KilometerPerSec;

define_unit!("Km/h", KilometerPerHour, Velocity, Kilometer::RATIO/Hour::RATIO);
pub type KilometersPerHour = KilometerPerHour;

define_unit!("au/day", AUPerDay, Velocity, Au::RATIO/Day::RATIO);
pub type AUsPerDay = Quantity<AUPerDay>;


impl std::ops::Div<Days> for AstronomicalUnits {
    type Output = AUsPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        Self::Output::new(self.0 / rhs.0)
    }
}

impl std::ops::Mul<Days> for AUsPerDay {
    type Output = AstronomicalUnits;

    fn mul(self, rhs: Days) -> Self::Output {
        Self::Output::new(self.0 * rhs.value())
    }
}
