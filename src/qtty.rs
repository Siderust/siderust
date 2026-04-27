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

pub use crate::ext_qtty::{
    accel, angular_rate, velocity, Acceleration, AmountOfSubstance, Angular,
    AngularRate, Current, Dimension, Dimensionless, Energy, Force, Length,
    LuminousIntensity, Mass, Per, Power, Prod, Quantity, Real, Scalar, Temperature, Time,
    Transcendental, Unit, Velocity, Volume,
};
pub use crate::ext_qtty::angular_rate::AngularRateUnit;

pub use crate::ext_qtty::unit;
pub use crate::ext_qtty::{
    acceleration, angular, area, energy, force, length, mass, power, time, volume,
};

pub use crate::ext_qtty::acceleration::*;
pub use crate::ext_qtty::angular::*;
pub use crate::ext_qtty::area::*;
pub use crate::ext_qtty::energy::*;
pub use crate::ext_qtty::force::*;
pub use crate::ext_qtty::length::*;
pub use crate::ext_qtty::mass::*;
pub use crate::ext_qtty::power::*;
pub use crate::ext_qtty::time::*;
pub use crate::ext_qtty::volume::*;
