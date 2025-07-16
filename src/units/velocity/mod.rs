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

//mod dms_per_year;
mod au_per_day;

//pub use dms_per_year::*;
pub use au_per_day::*;

use crate::units::*;
use std::f64::consts::PI;

pub enum Velocity {}
impl Dimension for Velocity {}
pub trait VelocityUnit: Unit<Dim = Velocity> {}
impl<T: Unit<Dim = Velocity>> VelocityUnit for T {}

define_unit!("deg/day", DegreePerDay, Velocity, 1.0);
pub type DegreesPerDay = Quantity<DegreePerDay>;


define_unit!("deg/year", DegreePerYear, Velocity, 1.0/Year::RATIO);
pub type DegreesPerYear = Quantity<DegreePerYear>;

define_unit!("rad/day", RadianPerDay, Velocity, 180.0/PI);
pub type RadiansPerDay = Quantity<RadianPerDay>;


impl std::ops::Div<Days> for Radians {
    type Output = RadiansPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        RadiansPerDay::new(self.value() / rhs.value())
    }
}

impl std::ops::Mul<Days> for RadiansPerDay {
    type Output = Radians;

    fn mul(self, rhs: Days) -> Self::Output {
        Radians::new(self.0 * rhs.value())
    }
}


impl std::ops::Div<Years> for Degrees {
    type Output = DegreesPerYear;

    fn div(self, rhs: Years) -> Self::Output {
        // 1 year = 365.25 days
        DegreesPerYear::new(self.value() / rhs.value())
    }
}

impl std::ops::Mul<Years> for DegreesPerYear {
    type Output = Degrees;

    fn mul(self, rhs: Years) -> Self::Output {
        Degrees::new(self.0 * rhs.value())
    }
}
