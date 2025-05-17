//! Solar Radius (R☉) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Solar Radii (R☉)
//! and conversions to and from Kilometers (km).

use super::Kilometers;

pub const SOLAR_RADIUS: SolarRadius = SolarRadius::new(1.0);

/// A strongly-typed representation of a length in Solar Radii (R☉).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SolarRadius(f64);

impl SolarRadius {
    /// Number of kilometers in one Solar Radius (R☉).
    pub const KM_PER_SOLAR_RADIUS: f64 = 695_700.0;

    /// Creates a new `SolarRadius` from a value in R☉.
    ///
    /// # Arguments
    /// * `value` - The length in solar radii.
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner value in R☉.
    pub const fn value(&self) -> f64 {
        self.0
    }

    /// Converts this value to kilometers.
    pub const fn to_km(&self) -> Kilometers {
        Kilometers::new(self.0 * Self::KM_PER_SOLAR_RADIUS)
    }
}

/// Converts a `SolarRadius` to `Kilometers`.
impl From<SolarRadius> for Kilometers {
    fn from(r: SolarRadius) -> Self {
        r.to_km()
    }
}

/// Converts a `Kilometers` to `SolarRadius`.
impl From<Kilometers> for SolarRadius {
    fn from(km: Kilometers) -> Self {
        Self::new(km.value() / SolarRadius::KM_PER_SOLAR_RADIUS)
    }
}

impl std::fmt::Display for SolarRadius {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} R☉", self.value())
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(SolarRadius);
