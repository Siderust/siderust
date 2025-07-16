//! # Angular Units Module
//!
//! This module defines and handles various angular measurement units
//! used in astronomical calculations. It provides types and utilities
//! for converting between common angular formats such as:
//!
//! - **Degrees**: Decimal representation of an Angular.
//! - **Radians**: Standard mathematical unit for Angular measurement.
//! - **DMS (Degrees, Minutes, Seconds)**: Sexagesimal format used in many
//!   astronomical and geodetic contexts.
//! - **HMS (Hours, Minutes, Seconds)**: Commonly used in right ascension
//!   and time-based angular measurements.
//! - **Arcseconds**: Subdivision of degrees, used for high-precision measurements.
//!
//! ## Features
//! - Conversion between units (e.g., degrees to radians, DMS to degrees, HMS to degrees).
//! - Safe and precise parsing and formatting utilities.
//! - Strongly typed structs for each angular unit to avoid accidental misuse.
//! - Arithmetic operations for compatible angular units.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::*;
//!
//! let deg = Degrees::new(180.0);
//! let rad = deg.to::<Radian>();
//! assert_eq!(rad.value(), std::f64::consts::PI);
//!
//! let dms = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
//! let hms = HMS::new(5, 30, 0.0);
//! ```
//!
//! This module aims to make astronomical Angular manipulations explicit,
//! correct, and ergonomic.

mod dms;
mod hms;

pub use dms::*;
pub use hms::*;

use crate::units::{define_unit, Quantity, Dimension, Unit};
use std::f64::consts::{TAU, PI};

pub enum Angular {}
impl Dimension for Angular {}
pub trait AngularUnit: Unit<Dim = Angular> {
    const FULL_TURN: f64;
    const HALF_TURN: f64;
    const QUARTED_TURN: f64;
}
impl<T: Unit<Dim = Angular>> AngularUnit for T {
    const FULL_TURN: f64 = Radians::new(TAU).to::<T>().value();
    const HALF_TURN: f64 = Radians::new(TAU).to::<T>().value() * 0.5;
    const QUARTED_TURN: f64 = Radians::new(TAU).to::<T>().value() * 0.25 ;
}

impl<U: AngularUnit + Copy> Quantity<U> {
    pub const TAU: Radians = Radians::new(TAU);

    #[inline]
    pub fn sin(&self) -> f64 {
        self.to::<Radian>().value().sin()
    }

    #[inline]
    pub fn cos(&self) -> f64 {
        self.to::<Radian>().value().cos()
    }

    #[inline]
    pub fn tan(&self) -> f64 {
        self.to::<Radian>().value().tan()
    }

    #[inline]
    pub fn sin_cos(&self) -> (f64, f64) {
        self.to::<Radian>().value().sin_cos()
    }

    /// Sign of the *raw numeric* in this unit (same semantics as f64::signum()).
    #[inline]
    pub const fn signum(self) -> f64 {
        self.value().signum()
    }

    #[inline]
    pub fn normalize(self) -> Self {
        self.wrap_pos()
    }

    #[inline]
    pub fn wrap_pos(self) -> Self {
        Self::new(self.value().rem_euclid(U::FULL_TURN))
    }

    /// Wrap to (-HALF_TURN, HALF_TURN].  (Your original `wrap_signed` semantics.)
    #[inline]
    pub fn wrap_signed(self) -> Self {
        let full = U::FULL_TURN;
        let half = 0.5 * full;
        let x = self.value();
        let y = (x + half).rem_euclid(full) - half;
        let norm = if y <= -half { y + full } else { y };
        Self::new(norm)
    }

    /// Wrap to [-HALF_TURN, HALF_TURN) (alternate boundary).
    #[inline]
    pub fn wrap_signed_lo(self) -> Self {
        let mut y = self.wrap_signed().value(); // now in (-half,half]
        let half = 0.5 * U::FULL_TURN;
        if y > half - 0.0 { // exact compare is fine; y==half means +half
            y -= U::FULL_TURN; // move +half to -half
        }
        Self::new(y)
    }

    /// "Latitude fold": map to [-QUARTER_TURN, +QUARTER_TURN].
    ///
    /// Equivalent to your `wrap_quarter_fold` for degrees.
    #[inline]
    pub fn wrap_quarter_fold(self) -> Self {
        let full = U::FULL_TURN;
        let half = 0.5 * full;
        let quarter = 0.25 * full;
        let y = (self.value() + quarter).rem_euclid(full);
        // quarter - |y - half| yields [-quarter, quarter]
        Self::new(quarter - (y - half).abs())
    }

    /// Signed smallest separation in (-HALF_TURN, HALF_TURN].
    #[inline]
    pub fn signed_separation(self, other: Self) -> Self {
        (self - other).wrap_signed() // requires Sub impl returning Self
    }

    //#[inline]
    pub fn abs_separation(self, other: Self) -> Self {
        let sep = self.signed_separation(other);
        Self::new(sep.value().abs())
    }

}

define_unit!("Deg", Degree, Angular, 1.0);
pub type Deg = Degree;
pub type Degrees = Quantity<Deg>;
pub const DEG: Degrees = Degrees::new(1.0);

define_unit!("Rad", Radian, Angular, 180.0/PI);
pub type Rad = Radian;
pub type Radians = Quantity<Rad>;
pub const RAD: Radians = Radians::new(1.0);

define_unit!("Arcs", Arcsecond, Angular, 1.0/3600.0);
pub type Arcs = Arcsecond;
pub type Arcseconds = Quantity<Arcs>;
pub const ARCS: Arcseconds = Arcseconds::new(1.0);


impl Degrees {
        /// Returns the inner `f64` value for this Angular in degrees.
    pub const fn from_hms(hours: i32, minutes: u32, seconds: f64) -> Degrees {
        let h_deg = hours as f64 * 15.0;
        let m_deg = minutes as f64 * 15.0 / 60.0;
        let s_deg = seconds * 15.0 / 3600.0;
        Self::new(h_deg + m_deg + s_deg)
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_full_turn() {
        assert_eq!(Radian::FULL_TURN, TAU);
        assert_eq!(Degree::FULL_TURN, 360.0);
        assert_eq!(Arcsecond::FULL_TURN, 1_296_000.0);
    }
}
