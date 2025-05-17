//! Astronomical Unit (AU) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Astronomical Units (AU)
//! and conversions to and from Light Years (LY).

use super::Kilometers;

pub const AU: AstronomicalUnit = AstronomicalUnit::new(1.0);

/// A strongly-typed representation of a length in Astronomical Units (AU).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct AstronomicalUnit(f64);

impl AstronomicalUnit {

    pub const KM_PER_AU: f64 = 149_597_870.7;

    /// Creates a new `AstronomicalUnit` from a value in AU.
    ///
    /// # Arguments
    /// * `value` - The length in astronomical units.
    ///
    /// # Example
    /// ```
    /// use siderust::units::AstronomicalUnit;
    /// let au = AstronomicalUnit::new(1.0);
    /// ```
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    pub const fn to_light_year(&self) -> super::LightYear {
        super::LightYear::new(self.0 / super::LightYear::AU_PER_LY)
    }

    pub const fn to_km(&self) -> Kilometers {
        Kilometers::new(self.0 * Self::KM_PER_AU)
    }

    /// Returns the inner value in AU.
    pub const fn value(&self) -> f64 {
        self.0
    }
}

/// Converts an `AstronomicalUnit` to a `LightYear`.
impl From<AstronomicalUnit> for super::LightYear {
    fn from(ly: AstronomicalUnit) -> Self {
        ly.to_light_year()
    }
}

impl std::fmt::Display for AstronomicalUnit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} AU", self.0)
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(AstronomicalUnit);
