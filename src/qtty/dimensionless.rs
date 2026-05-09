// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Dimensionless quantity compatibility exports for siderust.
//!
//! Generic atmospheric and photometric dimensionless units are defined in the
//! canonical `qtty` crate and re-exported here for the existing `siderust::qtty`
//! surface. The IAU CIP coordinate remains local because it is specific to the
//! precession-nutation model rather than a general physical quantity.

pub use crate::ext_qtty::dimensionless::{
    Airmass, Airmasses, Albedo, Albedos, IlluminationFraction, IlluminationFractions, OpticalDepth,
    OpticalDepths, Ratio, Ratios, Refractivities, Refractivity, Transmittance, Transmittances,
};

use crate::ext_qtty::{Dimensionless, Quantity, Unit};

/// Unit marker for IAU 2006/2000A Celestial Intermediate Pole (CIP) X or Y
/// coordinate.
///
/// CIP X and Y are tabulated dimensionless small-angle quantities in the SOFA
/// framework. They remain in `siderust` because their semantics are tied to the
/// astronomy transform model rather than to general quantity algebra.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CipCoordinate;

impl Unit for CipCoordinate {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// `Quantity<CipCoordinate>` — a typed CIP X/Y coordinate scalar.
pub type CipCoordinates = Quantity<CipCoordinate>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reexported_dimensionless_units_round_trip() {
        assert_eq!(OpticalDepths::new(0.5_f64).value(), 0.5_f64);
        assert_eq!(Airmasses::new(1.0_f64).value(), 1.0_f64);
        assert_eq!(Transmittances::new(0.7_f64).value(), 0.7_f64);
        assert_eq!(Albedos::new(0.3_f64).value(), 0.3_f64);
        assert_eq!(IlluminationFractions::new(0.75_f64).value(), 0.75_f64);
        assert_eq!(Refractivities::new(2.7e-4_f64).value(), 2.7e-4_f64);
    }

    #[test]
    fn cip_coordinate_round_trip() {
        let q = CipCoordinates::new(-1.234e-3_f64);
        assert_eq!(q.value(), -1.234e-3_f64);
    }
}
