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
//! use siderust::units::{AU, LightYear};
//!
//! let au = 1.0*AU;
//! //let ly = LightYear::from(au); TODO
//! //assert!((ly.value() - 1.582e-5).abs() < 1e-8);
//! ```

mod solar_radius;
pub use solar_radius::*;



use paste::paste;
use crate::units::Quantity;

macro_rules! define_unit {
    ($symbol:ident, $name:ident, $dim:ty, $ratio:expr) => {

        #[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
        pub enum $name {}
        impl Unit for $name {
            const RATIO: f64 = $ratio;
            type Dim = $dim;
            const SYMBOL: &'static str = stringify!($symbol);
        }
        paste! {
            pub type [<$name s>] = Quantity<$name>;
            pub type $symbol = [<$name s>];
            pub const $symbol: [<$name s>] = [<$name s>]::new(1.0);
        }

        impl std::fmt::Display for Quantity<$name> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.value(), <$name>::SYMBOL)
            }
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

/// AstronomicalUnit -> LightYear.
impl From<AU> for LY { fn from(au: AU) -> Self { au.to::<LightYear>() } }

/// LightYear -> AstronomicalUnits.
impl From<LY> for AU { fn from(ly: LY) -> Self { ly.to::<AstronomicalUnit>() } }

// TODO: Remove me
impl Unit for f64 {
    const RATIO: f64 = 1.0;
    type Dim = Length;
    const SYMBOL: &'static str = "";
}
