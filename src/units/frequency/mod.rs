//! # Frequency Units Module
//!
//! This module provides types and utilities for handling frequency-related calculations
//! in astronomical and scientific contexts. It includes representations for angular
//! frequency such as radians per day and conversions from angular and time units.
//!
//! ## Features
//! - **RadiansPerDay**: Represents angular frequency in radians per day.
//! - Arithmetic operations for frequency types.
//! - Conversion from `Radians` and `Days`.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{Radians, Days, RadiansPerDay};
//!
//! let angle = Radians::new(std::f64::consts::PI);
//! let time = Days::new(2.0);
//! let frequency = angle / time; // RadiansPerDay
//! assert_eq!(frequency.value(), std::f64::consts::PI / 2.0);
//! ```

use crate::units::*;
use std::f64::consts::PI;

pub enum Frequency {}
impl Dimension for Frequency {}
pub trait FrequencyUnit: Unit<Dim = Frequency> {}
impl<T: Unit<Dim = Frequency>> FrequencyUnit for T {}

define_unit!("deg/day", DegreePerDay, Frequency, 1.0);
pub type DegreesPerDay = Quantity<DegreePerDay>;


define_unit!("deg/year", DegreePerYear, Frequency, 1.0/Year::RATIO);
pub type DegreesPerYear = Quantity<DegreePerYear>;

define_unit!("rad/day", RadianPerDay, Frequency, 180.0/PI);
pub type RadiansPerDay = Quantity<RadianPerDay>;


impl std::ops::Div<Days> for Radians {
    type Output = RadiansPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        Self::Output::new(self.value() / rhs.value())
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
