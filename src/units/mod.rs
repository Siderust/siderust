//! # Units Module
//!
//! This module provides a comprehensive set of strongly-typed units and utilities
//! for astronomical and scientific calculations. It is designed to ensure correctness,
//! clarity, and ease of use when working with various units of measurement.
//!
//! ## Features
//! - **Time Units**: Includes representations for Days, Years, Julian Years, and Centuries.
//! - **Angular Units**: Provides types for Degrees, Radians, DMS (Degrees, Minutes, Seconds), HMS (HourAngles, Minutes, Seconds), and Arcseconds.
//! - **Length Units**: Includes types for meters and astronomical units (AstronomicalUnits).
//! - **Velocity Units**: Provides types for meters per second and kilometers per second.
//! - **Mass Units**: Includes types for kilograms and solar masses.
//! - **Power Units**: Includes types for watts and solar luminosity.
//! - **Arithmetic Operations**: Supports arithmetic operations between compatible units, ensuring type safety.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::*;
//!
//! // Angular Units
//! let degrees = Degrees::new(180.0);
//! let radians = degrees.to::<Radian>();
//! let dms = Degrees::from_dms(12, 34, 56.0);
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

pub mod angular;
pub mod frequency;
pub mod length;
pub mod mass;
pub mod power;
pub mod time;
pub mod unitless;
pub mod velocity;

pub use angular::*;
pub use frequency::*;
pub use length::*;
pub use mass::*;
pub use power::*;
pub use time::*;
pub use unitless::*;
pub use velocity::*;

use core::cmp::*;
use core::marker::PhantomData;
use core::ops::*;
use std::fmt::*;

/// Marker trait for **dimensions** (Length, Time, Mass …).
///
/// A *dimension* is the category that distinguishes a metre from a second.
/// You usually model each dimension as an empty enum:
///
/// ```rust
/// use siderust::units::Dimension;
/// #[derive(Debug)]
/// pub enum Length {}
/// impl Dimension for Length {}
/// ```
pub trait Dimension {}

/// Trait implemented by every **unit** type generated through [`define_unit!`].
///
/// * `RATIO` expresses how many of this unit fit into the *canonical* unit
///   of the same dimension.
///   Example: The ratio for centimetres is `1.000` because `1 km = 1.000 m`.
///
/// * `SYMBOL` is the printable string (e.g. `"m"` or `"cm"`).
///
/// * `Dim` ties the unit to its underlying [`Dimension`].
///
/// # Safety
/// The trait is `Copy + 'static`, so types must be zero-sized marker enums.
pub trait Unit: Copy + PartialEq + Debug + 'static {
    /// Unit-to-canonical conversion factor.
    const RATIO: f64;

    /// Dimension to which this unit belongs.
    type Dim: Dimension;

    /// Printable symbol, shown by [`Display`](std::fmt::Display).
    const SYMBOL: &'static str;
}

/// Dimension formed by dividing one [`Dimension`] by another.
///
/// This is used to model composite dimensions such as `Length/Time`
/// for velocities or `Angular/Time` for frequencies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DivDim<N: Dimension, D: Dimension>(PhantomData<(N, D)>);
impl<N: Dimension, D: Dimension> Dimension for DivDim<N, D> {}

/// Unit representing the division of two other units.
///
/// `Per<N, D>` corresponds to `N / D` and carries both the
/// dimensional information and the scaling ratio between the
/// constituent units. It is generic over any numerator and
/// denominator units, which allows implementing arithmetic
/// generically for all pairs without bespoke macros.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Per<N: Unit, D: Unit>(PhantomData<(N, D)>);

impl<N: Unit, D: Unit> Unit for Per<N, D> {
    const RATIO: f64 = N::RATIO / D::RATIO;
    type Dim = DivDim<N::Dim, D::Dim>;
    // The symbol is constructed at formatting time since generic
    // constants cannot concatenate at compile time.
    const SYMBOL: &'static str = "";
}

impl<N: Unit, D: Unit> fmt::Display for Quantity<Per<N, D>> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}/{}", self.value(), N::SYMBOL, D::SYMBOL)
    }
}

/// Numeric value tagged with a unit at compile time.
///
/// Arithmetic between mismatched dimensions is a **compile-time error**.
/// Scalar factors (`f64`) are allowed on either side of `*`/`/`.
///
/// ```rust
/// use siderust::units::length::Meter;
/// use siderust::units::time::Second;
/// use siderust::units::Quantity;
///
/// let speed = Quantity::<Meter>::new(10.0) / Quantity::<Second>::new(2.0);
/// //             ^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
/// //            10 m (Length)                2 s (Time)
///
/// // error[E0308]: binary operation `/` cannot be applied to two `Quantity<_>`
/// ```
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Quantity<U: Unit>(f64, PhantomData<U>);

impl<U: Unit + Copy> Quantity<U> {
    pub const NAN: Self = Self::new(f64::NAN);

    pub const fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    pub const fn value(self) -> f64 {
        self.0
    }

    pub fn abs(self) -> Self {
        Self::new(self.0.abs())
    }

    pub const fn to<T: Unit<Dim = U::Dim>>(self) -> Quantity<T> {
        Quantity::<T>::new(self.0 * (U::RATIO / T::RATIO))
    }

    pub const fn min(&self, other: Quantity<U>) -> Quantity<U> {
        Quantity::<U>::new(self.value().min(other.value()))
    }
}

impl<U> Add for Quantity<U>
where
    U: Unit,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.0 + rhs.0)
    }
}

impl<U> AddAssign for Quantity<U>
where
    U: Unit,
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
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.0 - rhs.0)
    }
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
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.0 * rhs)
    }
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
    fn div(self, rhs: f64) -> Self {
        Self::new(self.0 / rhs)
    }
}

impl<N: Unit, D: Unit> Div<Quantity<D>> for Quantity<N> {
    type Output = Quantity<Per<N, D>>;

    fn div(self, rhs: Quantity<D>) -> Self::Output {
        Quantity::<Per<N, D>>::new(self.0 / rhs.0)
    }
}

impl<N: Unit, D: Unit> Mul<Quantity<D>> for Quantity<Per<N, D>> {
    type Output = Quantity<N>;

    fn mul(self, rhs: Quantity<D>) -> Self::Output {
        Quantity::<N>::new(self.0 * rhs.value())
    }
}

impl<N: Unit, D: Unit> Mul<Quantity<Per<N, D>>> for Quantity<D> {
    type Output = Quantity<N>;

    fn mul(self, rhs: Quantity<Per<N, D>>) -> Self::Output {
        rhs * self
    }
}

impl<N: Unit, D: Unit> Div<Quantity<Per<N, D>>> for Quantity<N> {
    type Output = Quantity<D>;

    fn div(self, rhs: Quantity<Per<N, D>>) -> Self::Output {
        Quantity::<D>::new(self.0 / rhs.0)
    }
}

impl<U> DivAssign for Quantity<U>
where
    U: Unit,
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
    fn rem(self, rhs: f64) -> Self {
        Self::new(self.0 % rhs)
    }
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
    fn neg(self) -> Self {
        Self::new(-self.0)
    }
}

impl<U> From<f64> for Quantity<U>
where
    U: Unit,
{
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}
/// Generate a **unit type** and its [`Display`] implementation.
#[macro_export]
macro_rules! define_unit {
    ($symbol:expr, $name:ident, $dim:ty, $ratio:expr) => {
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

use define_unit;
