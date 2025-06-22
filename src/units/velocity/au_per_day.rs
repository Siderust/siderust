//! AUsPerDay type and conversions.
//!
//! Provides a strongly-typed representation of linear velocity in Astronomical Uits (AU) per day.

use crate::units::time::Days;
use crate::units::length::AstronomicalUnit;

/// A strongly-typed representation of linear velocity in AstronomicalUnit per day.
#[derive(Debug, Clone, PartialEq)]
pub struct AUsPerDay(pub AstronomicalUnit, pub Days);

impl AUsPerDay {
    /// Creates a new `AUsPerDay` from a value in AstronomicalUnit per day.
    pub fn new(aus: AstronomicalUnit, days: Days) -> Self {
        Self(aus, days)
    }

    pub fn from_decimal(aus_per_day: f64) -> Self {
        Self(AstronomicalUnit::new(aus_per_day), Days::new(1.0))
    }

    /// Returns the f64 result of AstronomicalUnit per day.
    pub fn value(&self) -> f64 {
        self.0.value() / self.1.value()
    }
}

impl std::fmt::Display for AUsPerDay {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.6} AU/day", self.value())
    }
}

/// Division of `AstronomicalUnit` by `Days` yields `AUsPerDay`.
impl std::ops::Div<Days> for AstronomicalUnit {
    type Output = AUsPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        AUsPerDay::new(self, rhs)
    }
}

/// Multiplication of `AUsPerDay` by `Days` yields `AstronomicalUnit`.
impl std::ops::Mul<Days> for AUsPerDay {
    type Output = AstronomicalUnit;

    fn mul(self, rhs: Days) -> Self::Output {
        self.0 * (rhs / self.1)
    }
}
