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

use super::*;
use paste::paste;

pub enum Velocity {}
impl Dimension for Velocity {}
pub trait VelocityUnit: Unit<Dim = Velocity> {}
impl<T: Unit<Dim = Velocity>> VelocityUnit for T {}

macro_rules! velocity_unit_auto {
    ($Num:ident, $Den:ident) => {
        paste! {
            define_unit!(
                concat!($Num::SYMBOL, "/", $Den::SYMBOL),
                [<$Num Per $Den>],
                Velocity,
                $Num::RATIO / $Den::RATIO
            );
            pub type [<$Num sPer $Den>] = Quantity<[<$Num Per $Den>]>;
            
            impl std::ops::Div<Quantity<$Den>> for Quantity<$Num> {
                type Output = Quantity<[<$Num Per $Den>]>;

                fn div(self, rhs: Quantity<$Den>) -> Self::Output {
                    Self::Output::new(self.0 / rhs.0)
                }
            }

            impl std::ops::Mul<Quantity<$Den>> for  Quantity<[<$Num Per $Den>]> {
                type Output = Quantity<$Num>;

                fn mul(self, rhs: Quantity<$Den>) -> Self::Output {
                    Self::Output::new(self.0 * rhs.value())
                }
            }

        }
    };
}

velocity_unit_auto!(Meter, Second);      // "m/s"
velocity_unit_auto!(Kilometer, Second);  // "Km/s"
velocity_unit_auto!(Kilometer, Hour);    // "Km/h"
velocity_unit_auto!(Au, Day);            // "au/day"

