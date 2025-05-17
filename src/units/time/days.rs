//! Time unit: Days
//!
//! This module provides the `Days` struct, representing a span of time in days.
//! It includes constants and arithmetic operations for working with days.

pub const DAY: Days = Days::new(1.0);

/// Represents a span of time in days.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Days(f64);


/// Number of seconds in a day (24h)
pub const SECONDS_PER_DAY: f64 = 86_400.0;

impl Days {
    pub const IN_A_YEAR: Days = Days(365.25);

    pub const fn new(value: f64) -> Self {
        Days(value)
    }

    #[inline]
    pub fn value(&self) -> f64 {
        self.0
    }
}

impl std::fmt::Display for Days {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} Days", self.0)
    }
}

impl std::ops::Div<Days> for Days {
    type Output = f64;
    fn div(self, days: Days) -> f64 {
        self.value() / days.value()
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(Days);
