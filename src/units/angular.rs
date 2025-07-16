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
//! let dms = Degrees::from_dms(12, 34, 56.0);
//! let hms = HourAngles::from_hms(5, 30, 0.0);
//! ```
//!
//! This module aims to make astronomical Angular manipulations explicit,
//! correct, and ergonomic.

use super::*;
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

define_unit!("Mas", MilliArcsecond, Angular, 1.0 / 3_600_000.0);
pub type Mas = MilliArcsecond;
pub type MilliArcseconds = Quantity<Mas>;
pub const MAS: MilliArcseconds = MilliArcseconds::new(1.0);


define_unit!("Hms", HourAngle, Angular, 15.0);
pub type Hms = HourAngle;
pub type HourAngles = Quantity<Hms>;
pub const HOUR_ANGLE: HourAngles = HourAngles::new(1.0);

impl HourAngles {
    pub const fn from_hms(hours: i32, minutes: u32, seconds: f64) -> Self {
        let sign = if hours < 0 { -1.0 } else { 1.0 };
        let h_abs = if hours < 0 { -hours } else { hours } as f64;
        let m = minutes as f64 / 60.0;
        let s = seconds / 3600.0;
        let total_hours = sign * (h_abs + m + s);
        Self::new(total_hours)
    }
}

impl Degrees {
    /// Construct from DMS components. Sign is taken from `deg`.
    /// `min` and `sec` are treated as magnitudes (>=0); no range check is enforced.
    /// Use `.normalize()` or `.wrap_*()` after constructing if you need canonical ranges.
    pub const fn from_dms(deg: i32, min: u32, sec: f64) -> Self {
        let sign = if deg < 0 { -1.0 } else { 1.0 };
        let d_abs = if deg < 0 { -deg } else { deg } as f64;
        let m = min as f64 / 60.0;
        let s = sec / 3600.0;
        let total = sign * (d_abs + m + s);
        Self::new(total)
    }

    /// Construct from explicit sign and magnitude components.
    /// `sign` should be -1, 0, or +1 (0 treated as +1 unless all components zero).
    pub const fn from_dms_sign(sign: i8, deg: u32, min: u32, sec: f64) -> Self {
        let s = if sign < 0 { -1.0 } else { 1.0 };
        let total = (deg as f64) + (min as f64)/60.0 + (sec/3600.0);
        Self::new(s * total)
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
