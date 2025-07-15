//! # Units Module
//!
//! This module provides a comprehensive set of strongly-typed units and utilities
//! for astronomical and scientific calculations. It is designed to ensure correctness,
//! clarity, and ease of use when working with various units of measurement.
//!
//! ## Features
//! - **Time Units**: Includes representations for Days, Years, Julian Years, and Centuries.
//! - **Angular Units**: Provides types for Degrees, Radians, DMS (Degrees, Minutes, Seconds), HMS (Hours, Minutes, Seconds), and Arcseconds.
//! - **Length Units**: Includes types for meters and astronomical units (AstronomicalUnits).
//! - **Velocity Units**: Provides types for meters per second and kilometers per second.
//! - **Mass Units**: Includes types for kilograms and solar masses.
//! - **Power Units**: Includes types for watts and solar luminosity.
//! - **Arithmetic Operations**: Supports arithmetic operations between compatible units, ensuring type safety.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{Days, Degrees, Radians, DMS, Kilograms, SolarMasses};
//!
//! // Angular Units
//! let degrees = Degrees::new(180.0);
//! let radians = Radians::from(degrees);
//! let dms = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
//!
//! // Mass Units
//! let mass_kg = Kilograms::new(5.0);
//! let mass_solar = SolarMasses::new(2.0);
//!
//! // Conversions
//! let dms_to_decimal = dms.value();
//!
//! assert_eq!(radians.value(), std::f64::consts::PI);
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

use core::marker::PhantomData;
use core::ops::*;
use core::cmp::*;
use std::fmt::*;

pub trait Dimension {} // (Length, Time, Mass…).

pub trait Unit: Copy + PartialEq + Debug + 'static {
    const RATIO: f64;
    type Dim: Dimension;
    const SYMBOL: &'static str;
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Quantity<U: Unit>(f64, PhantomData<U>);

impl<U: Unit + Copy> Quantity<U> {
    pub const NAN: Self = Self::new(f64::NAN);

    pub const fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    pub const fn value(self) -> f64 { self.0 }

    pub fn abs(self) -> Self { Self::new(self.0.abs()) }

    pub fn to<T: Unit<Dim = U::Dim>>(self) -> Quantity<T> {
        Quantity::<T>::new(self.0 * (U::RATIO / T::RATIO))
    }
}


impl<U> Add for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self { Self::new(self.0 + rhs.0) }
}

impl<U> AddAssign for Quantity<U>
where
    U: Unit
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<U> Sub for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { Self::new(self.0 - rhs.0) }
}

impl<U> SubAssign for Quantity<U>
where
    U: Unit,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl<U> Mul<f64> for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self { Self::new(self.0 * rhs) }
}

impl<U> Mul<Quantity<U>> for f64
where
    U: Unit,
{
    type Output = Quantity<U>;
    fn mul(self, rhs: Quantity<U>) -> Self::Output {
        rhs * self
    }
}

impl<U> Div<f64> for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn div(self, rhs: f64) -> Self { Self::new(self.0 / rhs) }
}

impl<U> Div<Quantity<U>> for Quantity<U>
where
    U: Unit,
{
    type Output = f64;
    fn div(self, rhs: Quantity<U>) -> Self::Output { self.0 / rhs.0 }
}

impl<U> DivAssign for Quantity<U>
where
    U: Unit
{
    fn div_assign(&mut self, rhs: Self) {
        self.0 /= rhs.0;
    }
}


impl<U> Rem<f64> for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn rem(self, rhs: f64) -> Self { Self::new(self.0 % rhs) }
}

impl<U> PartialEq<f64> for Quantity<U>
where
    U: Unit,
{
    fn eq(&self, other: &f64) -> bool {
        self.0 == *other
    }
}

impl<U> Neg for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.0) }
}

impl<U> From<f64> for Quantity<U>
where
    U: Unit,
{
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

#[macro_export]
macro_rules! define_unit {
    ($symbol:literal, $name:ident, $dim:ty, $ratio:expr) => {

        #[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
        pub enum $name {}
        impl Unit for $name {
            const RATIO: f64 = $ratio;
            type Dim = $dim;
            const SYMBOL: &'static str = stringify!($symbol);
        }
        impl std::fmt::Display for Quantity<$name> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.value(), <$name>::SYMBOL)
            }
        }
    };
}

pub(crate) use define_unit;