//! DmsPerYear type and conversions.
//!
//! Provides a strongly-typed representation of angular velocity in DMS (degrees, minutes, seconds) per year.

use crate::units::time::Years;
use crate::units::angular::{DMS, Degrees};

/// A strongly-typed representation of angular velocity in DMS per year.
#[derive(Debug, Clone, PartialEq)]
pub struct DmsPerYear(pub DMS, pub Years);

impl DmsPerYear {
    /// Creates a new `DmsPerYear` from a value in DMS per year.
    pub fn new(dms: DMS, years: Years) -> Self {
        Self(dms, years)
    }

    pub fn from_decimal(degrees_years: f64) -> Self {
        Self(DMS::from_degrees(Degrees::new(degrees_years)), Years::new(1.0))
    }

    /// Returns the f64 result of DMS per year.
    pub fn value(&self) -> f64 {
        self.0.value() / self.1.value()
    }
}

impl std::fmt::Display for DmsPerYear {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.6} DMS/year", self.value())
    }
}

/// Division of `DMS` by `Years` yields `DmsPerYear`.
impl std::ops::Div<Years> for DMS {
    type Output = DmsPerYear;

    fn div(self, rhs: Years) -> Self::Output {
        // 1 year = 365.25 days
        DmsPerYear::new(self, rhs)
    }
}

/// Multiplication of `DmsPerYear` by `Years` yields `DMS`.
impl std::ops::Mul<Years> for DmsPerYear {
    type Output = DMS;

    fn mul(self, rhs: Years) -> Self::Output {
        self.0 * (rhs / self.1)
    }
}
