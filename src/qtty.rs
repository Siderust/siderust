// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # `siderust::qtty` — typed physical quantities
//!
//! This module re-exports the external [`qtty`](https://crates.io/crates/qtty) crate
//! so you can write `siderust::qtty::Kilometer` without adding a separate
//! `qtty` dependency in your `Cargo.toml`. The pattern mirrors
//! `siderust::time` (a façade over `tempoch`).
//!
//! Adding `qtty` directly as a dependency gives you identical types; this
//! re-export is purely for ergonomics.
//!
//! ## Compatibility façade
//!
//! Older `siderust` code expects `use qtty::*;` to bring unit markers such as
//! `Day` and `AstronomicalUnit` into scope alongside plural quantity aliases
//! like `Days` and `AstronomicalUnits`. Modern `qtty` instead uses root
//! singular names for quantity aliases. Re-export the legacy mixed surface here
//! so the rest of this crate can keep its existing imports while depending on
//! newer `qtty`.

pub use crate::ext_qtty::angular_rate::AngularRateUnit;
pub use crate::ext_qtty::{
    accel, angular_rate, velocity, Acceleration, AmountOfSubstance, Angular, AngularRate, Current,
    Dimension, Dimensionless, Energy, Force, Length, LuminousIntensity, Mass, Per, Power, Prod,
    Quantity, Real, Scalar, Temperature, Time, Transcendental, Unit, Velocity, Volume,
};

pub use crate::ext_qtty::unit;
pub use crate::ext_qtty::{
    acceleration, angular, area, density, energy, force, length, mass, power, radiometry,
    solid_angle, time, volume,
};

pub use crate::ext_qtty::solid_angle::*;

pub use crate::ext_qtty::frequency;
pub use crate::ext_qtty::frequency::*;

pub use crate::ext_qtty::acceleration::*;
pub use crate::ext_qtty::angular::*;
pub use crate::ext_qtty::area::*;
pub use crate::ext_qtty::density::*;
pub use crate::ext_qtty::energy::*;
pub use crate::ext_qtty::force::*;
pub use crate::ext_qtty::length::*;
pub use crate::ext_qtty::mass::*;
pub use crate::ext_qtty::power::*;
pub use crate::ext_qtty::pressure;
pub use crate::ext_qtty::pressure::*;
pub use crate::ext_qtty::radiometry::*;
pub use crate::ext_qtty::temperature;
pub use crate::ext_qtty::temperature::*;
pub use crate::ext_qtty::time::*;
// Explicit import overrides the glob above: `crate::qtty::Second` = the
// quantity type `Quantity<unit::Second, f64>`, not the unit marker struct.
pub use crate::ext_qtty::volume::*;
pub use crate::ext_qtty::Second;

/// Re-export `GravitationalParameter` from `qtty::dynamics` so that
/// `crate::qtty::GravitationalParameter` is available in the `siderust` module tree.
pub use crate::ext_qtty::GravitationalParameter;

/// Re-export astrodynamics typed quantities used by `dynamics::state` and
/// force models.
pub use crate::ext_qtty::dynamics::{
    AreaToMass, AreaToMassUnit, DragCoefficient, InverseSecond, InverseSeconds, J2Coefficient,
    KmPerSecond, KmPerSecondSquared, KmPerSeconds, KmPerSecondsSquared, SrpCoefficient,
    SPEED_OF_LIGHT_KM_S,
};

/// Re-export standard gravitational parameter constants for common solar-system bodies.
///
/// All values are in km³/s² from JPL DE430 (or IAU 2012 for the Sun).
pub use crate::ext_qtty::dynamics::{GM_EARTH, GM_MOON, GM_SUN};

/// Standard gravitational parameter of Mercury (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_MERCURY: GravitationalParameter = GravitationalParameter::new(2.203_187_832_8e4);

/// Standard gravitational parameter of Venus (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_VENUS: GravitationalParameter = GravitationalParameter::new(3.248_585_920_0e5);

/// Standard gravitational parameter of Mars (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_MARS: GravitationalParameter = GravitationalParameter::new(4.282_837_581_6e4);

/// Standard gravitational parameter of Jupiter (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_JUPITER: GravitationalParameter = GravitationalParameter::new(1.267_127_648_0e8);

/// Standard gravitational parameter of Saturn (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_SATURN: GravitationalParameter = GravitationalParameter::new(3.794_058_520_0e7);

/// Standard gravitational parameter of Uranus (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_URANUS: GravitationalParameter = GravitationalParameter::new(5.794_548_600_0e6);

/// Standard gravitational parameter of Neptune (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_NEPTUNE: GravitationalParameter = GravitationalParameter::new(6.836_527_100_5e6);

/// Standard gravitational parameter of Pluto (system) in km³/s².
///
/// Source: JPL DE430.
pub const GM_PLUTO: GravitationalParameter = GravitationalParameter::new(9.770e2);

/// Re-export typed integrator tolerances for use in dynamics integrators.
pub use crate::ext_qtty::tolerances::IntegratorTolerances;

// ---------------------------------------------------------------------------
// Dimensionless quantities
// ---------------------------------------------------------------------------

// Generic atmospheric and photometric dimensionless units are re-exported from
// the canonical `qtty` crate. The IAU CIP coordinate stays local because its
// semantics are tied to the precession-nutation model, not general quantity algebra.

pub use crate::ext_qtty::dimensionless::{
    Airmass, Airmasses, Albedo, Albedos, IlluminationFraction, IlluminationFractions, OpticalDepth,
    OpticalDepths, Ratio, Ratios, Refractivities, Refractivity, Transmittance, Transmittances,
};

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
