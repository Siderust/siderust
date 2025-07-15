//! RadiansPerDay type and conversions.
//!
//! Provides a strongly-typed representation of angular velocity in radians per day,
//! and conversion from angular displacement and time.

use crate::units::{Days, Radians};

/// A strongly-typed representation of angular velocity in radians per day.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct RadiansPerDay(pub f64);

impl RadiansPerDay {
    /// Creates a new `RadiansPerDay` from a value in radians per day.
    ///
    /// # Arguments
    /// * `value` - The angular velocity in radians per day.
    ///
    /// # Example
    /// ```
    /// use siderust::units::RadiansPerDay;
    /// let v = RadiansPerDay::new(1.0);
    /// ```
    pub fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner value in radians per day.
    pub fn value(&self) -> f64 {
        self.0
    }
}

/// Implement `Display` for `RadiansPerDay`.
impl std::fmt::Display for RadiansPerDay {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.6} rad/day", self.0)
    }
}

/// Division of `Radians` by `Days` yields `RadiansPerDay`.
impl std::ops::Div<Days> for Radians {
    type Output = RadiansPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        RadiansPerDay::new(self.value() / rhs.value())
    }
}

/// Multiplication of `RadiansPerDay` per `Days` yields `Radians`.
impl std::ops::Mul<Days> for RadiansPerDay {
    type Output = Radians;

    fn mul(self, rhs: Days) -> Self::Output {
        Radians::new(self.0 * rhs.value())
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(RadiansPerDay);
