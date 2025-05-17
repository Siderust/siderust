//! Solar Mass (M☉) type and conversions.
//!
//! Provides a strongly-typed representation of a mass in Solar Masses (M☉).

pub const SOLAR_MASS: SolarMass = SolarMass::new(1.0);

/// A strongly-typed representation of a mass in Solar Masses (M☉).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SolarMass(f64);

impl SolarMass {
    pub const SOLAR_MASS_KG: f64 = 1.98847e30;

    /// Creates a new `SolarMass` from a value in solar masses.
    ///
    /// # Arguments
    /// * `value` - The mass in solar masses.
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner value in solar masses.
    pub const fn value(&self) -> f64 {
        self.0
    }
}

impl std::fmt::Display for SolarMass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} M☉", self.0)
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(SolarMass);
