//! AUsPerDay type and conversions.
//!
//! Provides a strongly-typed representation of linear velocity in Astronomical Uits (AstronomicalUnits) per day.

use crate::units::time::Days;
use crate::units::AstronomicalUnits;

/// A strongly-typed representation of linear velocity in AstronomicalUnits per day.
#[derive(Clone)]
pub struct AUsPerDay(pub AstronomicalUnits, pub Days);

impl AUsPerDay {
    /// Creates a new `AUsPerDay` from a value in AstronomicalUnits per day.
    pub fn new(aus: AstronomicalUnits, days: Days) -> Self {
        Self(aus, days)
    }

    pub fn from_decimal(aus_per_day: f64) -> Self {
        Self(AstronomicalUnits::new(aus_per_day), Days::new(1.0))
    }

    /// Returns the f64 result of AstronomicalUnits per day.
    pub fn value(&self) -> f64 {
        self.0.value() / self.1.value()
    }
}

impl std::fmt::Display for AUsPerDay {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.6} AstronomicalUnits/day", self.value())
    }
}

/// Division of `AstronomicalUnits` by `Days` yields `AUsPerDay`.
impl std::ops::Div<Days> for AstronomicalUnits {
    type Output = AUsPerDay;

    fn div(self, rhs: Days) -> Self::Output {
        AUsPerDay::new(self, rhs)
    }
}

/// Multiplication of `AUsPerDay` by `Days` yields `AstronomicalUnits`.
impl std::ops::Mul<Days> for AUsPerDay {
    type Output = AstronomicalUnits;

    fn mul(self, rhs: Days) -> Self::Output {
        self.0 * (rhs / self.1)
    }
}
