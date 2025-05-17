//! Solar Luminosity (L☉) type and conversions.
//!
//! Provides a strongly-typed representation of power in Solar Luminosities (L☉).

pub const SOLAR_LUMINOSITY: SolarLuminosity = SolarLuminosity::new(1.0);
/// A strongly-typed representation of power in Solar Luminosities (L☉).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SolarLuminosity(f64);

impl SolarLuminosity {
    /// The solar luminosity in SI units (watts).
    pub const SOLAR_LUMINOSITY_WATTS: f64 = 3.828e26;

    /// Creates a new `SolarLuminosity` from a value in solar luminosities.
    ///
    /// # Arguments
    /// * `value` - The power in solar luminosities.
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner value in solar luminosities.
    pub const fn value(&self) -> f64 {
        self.0
    }

    /// Converts this value to watts.
    pub fn to_watts(&self) -> f64 {
        self.0 * Self::SOLAR_LUMINOSITY_WATTS
    }

    /// Creates a `SolarLuminosity` from a value in watts.
    pub fn from_watts(watts: f64) -> Self {
        Self(watts / Self::SOLAR_LUMINOSITY_WATTS)
    }
}

impl std::fmt::Display for SolarLuminosity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} L☉", self.0)
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(SolarLuminosity);
