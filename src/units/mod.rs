//! Strongly typed physical quantities.
//!
//! The `units` module supplies zero-cost wrappers around primitive values to
//! encode their physical dimension in the type system. Arithmetic is only
//! permitted between compatible dimensions, reducing the risk of unit mistakes in
//! scientific code.
//!
//! ## Provided dimensions
//! - [`time`]: days, Julian years and related intervals.
//! - [`angular`]: degrees, radians, hour angles and arcseconds.
//! - [`length`]: metres, kilometres and astronomical units.
//! - [`velocity`], [`mass`], [`power`], [`frequency`] and [`unitless`].
//!
//! ## Example
//! ```rust
//! use siderust::units::*;
//! let angle = Degrees::new(180.0);
//! let radians = angle.to::<Radians>();
//! assert_eq!(radians.value(), std::f64::consts::PI);
//! ```

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

impl<N: Unit, D: Unit> Display for Quantity<Per<N, D>> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{} {}/{}", self.value(), N::SYMBOL, D::SYMBOL)
    }
}

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

impl<U: Unit> Add for Quantity<U> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.0 + rhs.0)
    }
}

impl<U: Unit> AddAssign for Quantity<U> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<U: Unit> Sub for Quantity<U> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.0 - rhs.0)
    }
}

impl<U: Unit> SubAssign for Quantity<U> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl<U: Unit> Mul<f64> for Quantity<U> {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.0 * rhs)
    }
}

impl<U: Unit> Mul<Quantity<U>> for f64 {
    type Output = Quantity<U>;
    fn mul(self, rhs: Quantity<U>) -> Self::Output {
        rhs * self
    }
}

impl<U: Unit> Div<f64> for Quantity<U> {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.0 / rhs)
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

impl<U: Unit> DivAssign for Quantity<U> {
    fn div_assign(&mut self, rhs: Self) {
        self.0 /= rhs.0;
    }
}

impl<U: Unit> Rem<f64> for Quantity<U> {
    type Output = Self;
    fn rem(self, rhs: f64) -> Self {
        Self::new(self.0 % rhs)
    }
}

impl<U: Unit> PartialEq<f64> for Quantity<U> {
    fn eq(&self, other: &f64) -> bool {
        self.0 == *other
    }
}

impl<U: Unit> Neg for Quantity<U> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.0)
    }
}

impl<U: Unit> From<f64> for Quantity<U> {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}


/* TODO: Requires specialization (nightly) see #16

impl<N: Unit, D: Unit> Div<Quantity<D>> for Quantity<N> {
    type Output = Quantity<Per<N, D>>;

    fn div(self, rhs: Quantity<D>) -> Self::Output {
        Quantity::<Per<N, D>>::new(self.0 / rhs.0)
    }
}

impl<N: Unit, D: Unit> Div<Quantity<Per<N, D>>> for Quantity<N> {
    type Output = Quantity<D>;

    fn div(self, rhs: Quantity<Per<N, D>>) -> Self::Output {
        Quantity::<D>::new(self.0 / rhs.0)
    }
}

*/

impl<N: Unit, D: Unit> std::ops::Div<Quantity<D>> for Quantity<N> {
    type Output = Quantity<Per<N, D>>;
    fn div(self, rhs: Quantity<D>) -> Self::Output {
        Quantity::new(self.value() / rhs.value())
    }
}

pub trait Simplify {
    type Out: Unit;
    fn simplify(self) -> Quantity<Self::Out>;
}

// U/U → Unitless
impl<U: Unit> Simplify for Quantity<Per<U, U>> {
    type Out = Unitless;
    fn simplify(self) -> Quantity<Unitless> {
        Quantity::new(self.value())
    }
}

// N / (N/D) → D
impl<N: Unit, D: Unit> Simplify for Quantity<Per<N, Per<N, D>>> {
    type Out = D;
    fn simplify(self) -> Quantity<D> {
        Quantity::new(self.value())
    }
}


impl<U: Unit> Quantity<Per<U, U>> {
    #[inline]
    pub fn asin(&self) -> f64 {
        self.value().asin()
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
