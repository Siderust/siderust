//! # Length Units Module
//!
//! This module provides types and utilities for handling length-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! length units and conversions between them.
//!
//! ## Features
//! - **Astronomical Unit (AU)**: The mean distance between the Earth and the Sun.
//! - **Light Year (LY)**: The distance light travels in one Julian year in vacuum.
//! - Conversion between AU and LY.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{AstronomicalUnit, LightYear};
//!
//! let au = AU::new(1.0);
//! let ly = LightYear::from(au);
//! assert!((ly.value() - 1.582e-5).abs() < 1e-8);
//! ```

mod astronomical_unit;
mod kilometer;
mod light_year;
mod solar_radius;

//pub use astronomical_unit::*;
//pub use kilometer::*;
//pub use light_year::*;
pub use solar_radius::*;



use paste::paste;
use crate::units::Quantity;

macro_rules! define_unit {
    ($abreviation:ident, $name:ident, $dim:ty, $ratio:expr) => {

        #[derive(Clone, Copy, Debug, PartialEq)]
        pub enum $name {}
        impl Unit for $name {
            const RATIO: f64 = $ratio;
            type Dim = $dim;
        }
        paste! {
            pub type [<$name s>] = Quantity<$name>;
            pub type $abreviation = [<$name s>];
            pub const $abreviation: [<$name s>] = [<$name s>]::new(1.0);
        }
    };
}

use crate::units::{Dimension, Unit};
pub enum Length {}
impl Dimension for Length {}
pub trait LengthUnit: Unit<Dim = Length> {}
impl<T: Unit<Dim = Length>> LengthUnit for T {}


define_unit!(M, Meter, Length, 1.0);
define_unit!(KM, Kilometer, Length, 1_000.0);
define_unit!(AU, AstronomicalUnit, Length, 149_597_870.7);
define_unit!(LY, LightYear, Length, 9_460_730_472_580.8);


// TODO: Remove me
impl Unit for f64 {
    const RATIO: f64 = 1.0;
    type Dim = Length;
}
