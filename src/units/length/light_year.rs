//! Light Year (LY) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Light Years (LY)
//! and conversions to and from Astronomical Units (AU).

use super::*;

impl LightYears {
    pub const AU_PER_LY: f64 = 63_241.1;
    pub const KM_PER_LY: f64 = 9_460_730_472_580.8;

    pub const fn to_km(&self) -> Kilometers {
        Kilometers::new(self.0 * Self::KM_PER_LY)
    }

    pub const fn to_au(&self) -> AstronomicalUnits {
        AU::new(self.0 * Self::AU_PER_LY)
    }

}

/// Converts a `LightYear` to `AstronomicalUnits`.
impl From<LightYears> for super::AstronomicalUnits {
    fn from(ly: LightYears) -> Self {
        ly.to_au()
    }
}
