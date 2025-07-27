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

use super::*;

pub enum Frequency {}
impl Dimension for Frequency {}
pub trait FrequencyUnit: Unit<Dim = Frequency> {}
impl<T: Unit<Dim = Frequency>> FrequencyUnit for T {}

macro_rules! frequency_unit_auto {
    ($Num:ident, $Den:ident) => {
        paste::paste! {
            define_unit!(
                concat!($Num::SYMBOL, "/", $Den::SYMBOL),
                [<$Num Per $Den>],
                Frequency,
                $Num::RATIO / $Den::RATIO
            );
            pub type [<$Num sPer $Den>] = Quantity<[<$Num Per $Den>]>;
            
            impl std::ops::Div<Quantity<$Den>> for Quantity<$Num> {
                type Output = Quantity<[<$Num Per $Den>]>;

                fn div(self, rhs: Quantity<$Den>) -> Self::Output {
                    Self::Output::new(self.0 / rhs.0)
                }
            }

            impl std::ops::Div<Quantity<[<$Num Per $Den>]>> for Quantity<$Num> {
                type Output = Quantity<$Den>;

                fn div(self, rhs: Quantity<[<$Num Per $Den>]>) -> Self::Output {
                    Self::Output::new(self.value() / rhs.value())
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

frequency_unit_auto!(Radian, Day);      // "rad/day"
frequency_unit_auto!(Degree, Day);      // "deg/day"
frequency_unit_auto!(Degree, Year);     // "deg/year"
frequency_unit_auto!(MilliArcsecond, Day); // "mas/day"
