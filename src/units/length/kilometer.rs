//! Kilometers (km) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Kilometers (km)
//! and conversions to and from Light Years (LY).

use super::{AstronomicalUnit, LightYear};

pub const KM: Kilometers = Kilometers::new(1.0);

/// A strongly-typed representation of a length in Kilometers (km).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Kilometers(f64);

impl Kilometers {

    /// Creates a new `Kilometers` from a value in km.
    ///
    /// # Arguments
    /// * `value` - The length in kilometers.
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    pub const fn to_au(&self) -> AstronomicalUnit {
        AstronomicalUnit::new(self.0 / AstronomicalUnit::KM_PER_AU)
    }

    /// Returns the inner value in kilometers.
    pub const fn value(&self) -> f64 {
        self.0
    }
}


/// Converts a `Kilometers` to `LightYear`.
impl From<Kilometers> for LightYear {
    fn from(km: Kilometers) -> Self {
        Self::new(km.value() / LightYear::KM_PER_LY)
    }
}

/// Converts a `Kilometers` to `AstronomicalUnits`.
impl From<Kilometers> for AstronomicalUnit {
    fn from(km: Kilometers) -> Self {
        km.to_au()
    }
}

impl std::fmt::Display for Kilometers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} km", self.value())
    }
}


// Import your macro path correctly if needed
crate::units::arithmetic_ops::impl_arithmetic_ops!(Kilometers);
