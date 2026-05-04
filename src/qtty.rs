// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Internal compatibility facade over the external `qtty` crate.
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
    acceleration, angular, area, energy, force, length, mass, power, radiometry, solid_angle, time,
    volume,
};

pub use crate::ext_qtty::solid_angle::*;

pub use crate::ext_qtty::acceleration::*;
pub use crate::ext_qtty::angular::*;
pub use crate::ext_qtty::area::*;
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
pub use crate::ext_qtty::Second;
pub use crate::ext_qtty::volume::*;
