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

use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Sub, SubAssign, Mul, Div, Rem, Neg};
use core::cmp::{PartialOrd, PartialEq};

pub trait Dimension {} // (Length, Time, Mass…).

pub trait Unit: Copy + 'static {
    const RATIO: f64;
    type Dim: Dimension;
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
        let base = self.0 * U::RATIO;
        Quantity::<T>::new(base / T::RATIO)
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


/*
pub trait Unit:
    Copy
    + Clone
    + PartialEq
    + PartialOrd
    + std::fmt::Debug
    + std::fmt::Display
    + std::ops::Add<Self, Output = Self>
    + std::ops::Sub<Self, Output = Self>
    + std::ops::Mul<f64, Output = Self>
    + std::ops::Div<f64, Output = Self>
    + From<f64>
    + Into<f64>
    + PartialEq<f64>
    + std::ops::AddAssign
    + std::ops::SubAssign
    + simba::scalar::ClosedAdd
    + simba::scalar::ClosedSub
    + nalgebra::Scalar
{
    const NAN: Self;

    fn sqrt(self) -> Self;
    fn powi(self, n: i32) -> Self;
    fn abs(self) -> Self;
}
// TODO: f64 has no Unit!! Remove me
impl Unit for f64 {
    const NAN: Self = f64::NAN;
    fn sqrt(self) -> Self { f64::sqrt(self) }
    fn powi(self, n: i32) -> Self { f64::powi(self, n) }
    fn abs(self) -> Self { f64::abs(self) }
}

macro_rules! impl_simple_unit {
    ($U:ty) => {
        impl Unit for $U {
            const NAN: Self = Self(f64::NAN);
            fn sqrt(self) -> Self { Self(self.0.sqrt()) }
            fn powi(self, n: i32) -> Self { Self(self.0.powi(n)) }
            fn abs(self) -> Self { Self(self.0.abs()) }
        }
    };
}

pub(crate) use impl_simple_unit;
*/