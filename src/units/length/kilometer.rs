//! Kilometers (km) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Kilometers (km)
//! and conversions to and from Light Years (LY).

use crate::units::*;

impl Kilometers {
    pub const fn to_au(&self) -> AstronomicalUnits {
        AU::new(self.0 / AstronomicalUnits::KM_PER_AU)
    }
}

