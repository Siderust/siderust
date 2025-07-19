//! Angular quantities and utilities.
//!
//! This module defines the **`Angular` dimension**, a blanket [`AngularUnit`] trait that extends
//! [`Unit`] for all angular units, common angular units (degrees, radians, arcseconds, etc.), and a set of
//! convenience methods on [`Quantity<U>`] where `U: AngularUnit`.
//!
//! # Design overview
//!
//! * **Canonical unit:** Degrees are taken as the canonical *scaling* unit for this dimension. That is,
//!   `Degree::RATIO == 1.0`, and all other angular units express how many *degrees* correspond to one of that unit.
//!   For example, `Radian::RATIO == 180.0 / PI` because 1 radian = 180/π degrees.
//! * **Associated constants:** The `AngularUnit` trait exposes precomputed constants (`FULL_TURN`, `HALF_TURN`,
//!   `QUARTED_TURN`) expressed *in the receiving unit* for ergonomic range‑wrapping.
//! * **Trigonometry:** `sin`, `cos`, `tan`, and `sin_cos` methods are provided on angular quantities; they convert to
//!   radians internally and then call the corresponding `f64` intrinsic.
//! * **Wrapping helpers:** Utility methods to wrap any angle into common ranges — `[0, 360)` (or unit equivalent),
//!   `(-180, 180]`, and the latitude‑style quarter fold `[-90, 90]`.
//!
//! ## Examples
//!
//! Convert between degrees and radians and evaluate a trig function:
//!
//! ```rust
//! use siderust::units::angular::{Degrees, Radians};
//! use siderust::units::Quantity;
//!
//! let angle: Degrees = Degrees::new(90.0);
//! let r: Radians = angle.to();
//! assert!((r.value() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
//! assert!((angle.sin() - 1.0).abs() < 1e-12);
//! ```
//!
//! Wrap into the conventional signed range:
//!
//! ```rust
//! use siderust::units::angular::Degrees;
//! let a = Degrees::new(370.0).wrap_signed();
//! assert_eq!(a.value(), 10.0);
//! ```

use super::*;
use std::f64::consts::{TAU, PI};

/// Dimension tag for angular measures (e.g., degrees, radians, arcseconds).
pub enum Angular {}
impl Dimension for Angular {}

/// Blanket extension trait for any [`Unit`] whose dimension is [`Angular`].
///
/// These associated constants provide the size of key turn fractions *expressed in the implementing unit*.
/// They are computed using a compiletime conversion from `TAU` radians (i.e., a full revolution).
///
/// > **Naming note:** The historical spelling `QUARTED_TURN` is retained for backward compatibility. It represents a
/// > quarter turn (90°).
pub trait AngularUnit: Unit<Dim = Angular> {
    const FULL_TURN: f64;
    const HALF_TURN: f64;
    const QUARTED_TURN: f64;
}
impl<T: Unit<Dim = Angular>> AngularUnit for T {
    /// One full revolution (360°) expressed in T unit.
    const FULL_TURN: f64 = Radians::new(TAU).to::<T>().value();
    /// Half a revolution (180°) expressed in T unit.
    const HALF_TURN: f64 = Radians::new(TAU).to::<T>().value() * 0.5;
    /// Quarter revolution (90°) expressed in T unit.
    const QUARTED_TURN: f64 = Radians::new(TAU).to::<T>().value() * 0.25 ;
}


impl<U: AngularUnit + Copy> Quantity<U> {
    /// Constant representing τ radians (2π rad == 360°).
    pub const TAU: Quantity::<U> = Quantity::<U>::new(U::FULL_TURN);

    /// Sine of the angle.
    #[inline]
    pub fn sin(&self) -> f64 {
        self.to::<Radian>().value().sin()
    }

    /// Cosine of the angle.
    #[inline]
    pub fn cos(&self) -> f64 {
        self.to::<Radian>().value().cos()
    }

    /// Tangent of the angle.
    #[inline]
    pub fn tan(&self) -> f64 {
        self.to::<Radian>().value().tan()
    }

    /// Simultaneously compute sine and cosine.
    #[inline]
    pub fn sin_cos(&self) -> (f64, f64) {
        self.to::<Radian>().value().sin_cos()
    }

    /// Sign of the *raw numeric* in this unit (same semantics as `f64::signum()`).
    #[inline]
    pub const fn signum(self) -> f64 {
        self.value().signum()
    }

    /// Normalize into the canonical positive range `[0, FULL_TURN)`.
    ///
    /// Shorthand for [`Self::wrap_pos`].
    #[inline]
    pub fn normalize(self) -> Self {
        self.wrap_pos()
    }

    /// Wrap into the positive range `[0, FULL_TURN)` using Euclidean remainder.
    #[inline]
    pub fn wrap_pos(self) -> Self {
        Self::new(self.value().rem_euclid(U::FULL_TURN))
    }

    /// Wrap into the signed range `(-HALF_TURN, HALF_TURN]`.
    ///
    /// *Upper bound is inclusive*; lower bound is exclusive. Useful for computing minimal signed angular differences.
    #[inline]
    pub fn wrap_signed(self) -> Self {
        let full = U::FULL_TURN;
        let half = 0.5 * full;
        let x = self.value();
        let y = (x + half).rem_euclid(full) - half;
        let norm = if y <= -half { y + full } else { y };
        Self::new(norm)
    }

    /// Wrap into the alternate signed range `[-HALF_TURN, HALF_TURN)`.
    ///
    /// Lower bound inclusive; upper bound exclusive. Equivalent to `self.wrap_signed()` with the boundary flipped.
    #[inline]
    pub fn wrap_signed_lo(self) -> Self {
        let mut y = self.wrap_signed().value(); // now in (-half, half]
        let half = 0.5 * U::FULL_TURN;
        if y > half { // move +half to -half
            y -= U::FULL_TURN;
        }
        Self::new(y)
    }

    /// "Latitude fold": map into `[-QUARTER_TURN, +QUARTER_TURN]`.
    ///
    /// Useful for folding polar coordinates (e.g., converting declination‑like angles to a limited range).
    #[inline]
    pub fn wrap_quarter_fold(self) -> Self {
        let full = U::FULL_TURN;
        let half = 0.5 * full;
        let quarter = 0.25 * full;
        let y = (self.value() + quarter).rem_euclid(full);
        // quarter - |y - half| yields [-quarter, quarter]
        Self::new(quarter - (y - half).abs())
    }

    /// Signed smallest angular separation in `(-HALF_TURN, HALF_TURN]`.
    #[inline]
    pub fn signed_separation(self, other: Self) -> Self {
        (self - other).wrap_signed()
    }

    /// Absolute smallest angular separation (magnitude only).
    #[inline]
    pub fn abs_separation(self, other: Self) -> Self {
        let sep = self.signed_separation(other);
        Self::new(sep.value().abs())
    }
}

define_unit!("Deg", Degree, Angular, 1.0);
/// Type alias shorthand for [`Degree`].
pub type Deg = Degree;
/// Convenience alias for a degree quantity.
pub type Degrees = Quantity<Deg>;
/// One degree.
pub const DEG: Degrees = Degrees::new(1.0);

// NOTE: 1 rad = 180/π degrees.
define_unit!("Rad", Radian, Angular, 180.0 / PI);
/// Type alias shorthand for [`Radian`].
pub type Rad = Radian;
/// Convenience alias for a radian quantity.
pub type Radians = Quantity<Rad>;
/// One radian.
pub const RAD: Radians = Radians::new(1.0);

// NOTE: 1 arcsecond = 1/3600 degree.
define_unit!("Arcs", Arcsecond, Angular, 1.0 / 3600.0);
/// Type alias shorthand for [`Arcsecond`].
pub type Arcs = Arcsecond;
/// Convenience alias for an arcsecond quantity.
pub type Arcseconds = Quantity<Arcs>;
/// One arcsecond.
pub const ARCS: Arcseconds = Arcseconds::new(1.0);

// NOTE: 1 milliarcsecond = 1/3_600_000 degree.
define_unit!("Mas", MilliArcsecond, Angular, 1.0 / 3_600_000.0);
/// Type alias shorthand for [`MilliArcsecond`].
pub type Mas = MilliArcsecond;
/// Convenience alias for a milliarcsecond quantity.
pub type MilliArcseconds = Quantity<Mas>;
/// One milliarcsecond.
pub const MAS: MilliArcseconds = MilliArcseconds::new(1.0);

// NOTE: 1 hour angle = 15 degrees.
define_unit!("Hms", HourAngle, Angular, 15.0);
/// Type alias shorthand for [`HourAngle`].
pub type Hms = HourAngle;
/// Convenience alias for an hour-angle quantity.
pub type HourAngles = Quantity<Hms>;
/// One hour angle hour (==15°).
pub const HOUR_ANGLE: HourAngles = HourAngles::new(1.0);

impl HourAngles {
    /// Construct from **HMS** components (`hours`, `minutes`, `seconds`).
    ///
    /// Sign is taken from `hours`; the `minutes` and `seconds` parameters are treated as magnitudes.
    ///
    /// ```rust
    /// use siderust::units::angular::HourAngles;
    /// let ra = HourAngles::from_hms(5, 30, 0.0); // 5h30m == 5.5h
    /// assert_eq!(ra.value(), 5.5);
    /// ```
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
    /// Construct from **DMS** components (`deg`, `min`, `sec`).
    ///
    /// Sign is taken from `deg`; the magnitude of `min` and `sec` is always added.
    /// No range checking is performed. Use one of the wrapping helpers if you need a canonical range.
    ///
    /// ```rust
    /// use siderust::units::angular::Degrees;
    /// let lat = Degrees::from_dms(-33, 52, 0.0); // −33°52′00″
    /// assert!(lat.value() < 0.0);
    /// ```
    pub const fn from_dms(deg: i32, min: u32, sec: f64) -> Self {
        let sign = if deg < 0 { -1.0 } else { 1.0 };
        let d_abs = if deg < 0 { -deg } else { deg } as f64;
        let m = min as f64 / 60.0;
        let s = sec / 3600.0;
        let total = sign * (d_abs + m + s);
        Self::new(total)
    }

    /// Construct from explicit sign and magnitude components.
    ///
    /// `sign` should be −1, 0, or +1 (0 treated as +1 unless all components are zero).
    pub const fn from_dms_sign(sign: i8, deg: u32, min: u32, sec: f64) -> Self {
        let s = if sign < 0 { -1.0 } else { 1.0 };
        let total = (deg as f64) + (min as f64)/60.0 + (sec/3600.0);
        Self::new(s * total)
    }
}

impl From<Degrees> for Radians { fn from(deg: Degrees) -> Self { deg.to::<Radian>() } }
impl From<Radians> for Degrees { fn from(rad: Radians) -> Self { rad.to::<Degree>() } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_full_turn() {
        assert_eq!(Radian::FULL_TURN, TAU);
        assert_eq!(Degree::FULL_TURN, 360.0);
        assert_eq!(Arcsecond::FULL_TURN, 1_296_000.0);
    }

    #[test]
    fn test_trig() {
        let a = Degrees::new(90.0);
        assert!((a.sin() - 1.0).abs() < 1e-12);
        assert!(a.cos().abs() < 1e-12);
    }

    #[test]
    fn test_wrap_signed() {
        let a = Degrees::new(370.0).wrap_signed();
        assert_eq!(a.value(), 10.0);
        let b = Degrees::new(-190.0).wrap_signed();
        assert_eq!(b.value(), 170.0);
    }
}
