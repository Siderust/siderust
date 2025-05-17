/// Number of seconds in a Julian year (365.25 days)
pub const SECONDS_PER_JULIAN_YEAR: f64 = 31_557_600.0;

/// Represents a time duration in Julian years.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct JulianYear(pub f64);

impl JulianYear {
    /// Create a new JulianYear from a value.
    pub fn new(years: f64) -> Self {
        JulianYear(years)
    }

    /// Get the value in Julian years.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Convert Julian years to seconds.
    pub fn to_seconds(&self) -> f64 {
        self.0 * SECONDS_PER_JULIAN_YEAR
    }

    /// Create a JulianYear from seconds.
    pub fn from_seconds(seconds: f64) -> Self {
        JulianYear(seconds / SECONDS_PER_JULIAN_YEAR)
    }
}

// Optionally, implement Display for pretty printing
use std::fmt;
impl fmt::Display for JulianYear {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} a (Julian)", self.0)
    }
}
