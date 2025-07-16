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
use std::f64::consts::PI;

pub enum Angular {}
impl Dimension for Angular {}
pub trait AngularUnit: Unit<Dim = Angular> {}
impl<T: Unit<Dim = Angular>> AngularUnit for T {}

impl<U: AngularUnit + Copy> Quantity<U> {
    /// Compute the sine of the Angular (in degrees), by converting internally to radians.
    #[inline]
    pub fn sin(&self) -> f64 {
        self.to::<Radian>().value().sin()
    }

    /// Compute the cosine of the Angular (in degrees).
    #[inline]
    pub fn cos(&self) -> f64 {
        self.to::<Radian>().value().cos()
    }

    /// Compute the tangent of the Angular (in degrees).
    #[inline]
    pub fn tan(&self) -> f64 {
        self.to::<Radian>().value().tan()
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

    /// Normalize an Angular in degrees to the range [0, 360).
    #[inline]
    pub fn normalize(self) -> Self {
        Self::new(self.value().rem_euclid(360.0))
    }

    /// Normalize to the range [−90, +90] in a symmetrical fashion.
    ///
    /// (Implementation depends on your use case—this is just an example.)
    #[inline]
    pub fn normalize_to_90_range(self) -> Self {
        let y = (self.value() + 90.0).rem_euclid(360.0);
        Self::new(90.0 - (y - 180.0).abs())
    }
    
    /// Normalize to the range [−180, +180].
    #[inline]
    pub fn normalize_to_180_range(self) -> Self {
        Self::new((self.value() + 180.0).rem_euclid(360.0) - 180.0)
    }

    pub fn diff_deg(self, other: Degrees) -> Self {
        // diferencia normalizada al intervalo 0 … 360
        let d = (self.value() - other.value()).rem_euclid(360.0);
        if d > 180.0 { Self::new(360.0 - d) } else { Self::new(d) }
    }
}


impl Radians {
    pub const TAU: Radians = Radians::new(std::f64::consts::TAU);

    /// Simultaneously computes the sine and cosine of the Angular.
    #[inline]
    pub fn sin_cos(self) -> (f64, f64) {
        self.value().sin_cos()
    }

    /// Returns the signum of the Angular (i.e. sign of the inner value).
    #[inline]
    pub fn signum(self) -> f64 {
        self.value().signum()
    }

    /// Normalize range (-π, π].
    #[inline]
    pub fn wrap_pi(self) -> Self {
        let x = self.value();
        let y = (x + PI).rem_euclid(2.0 * PI) - PI;
        let norm = if y <= -PI { y + 2.0 * PI } else { y };
        Radians::new(norm)
    }
}