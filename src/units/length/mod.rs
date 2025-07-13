//! # Length Units Module
//!
//! This module provides types and utilities for handling length-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! length units and conversions between them.
//!
//! ## Features
//! - **Astronomical Unit (AstronomicalUnits)**: The mean distance between the Earth and the Sun.
//! - **Light Year (LightYears)**: The distance light travels in one Julian year in vacuum.
//! - Conversion between AstronomicalUnits and LightYears.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{AU, LightYears};
//!
//! let au = 1.0*AU;
//! //let ly = LightYears::from(au); TODO
//! //assert!((ly.value() - 1.582e-5).abs() < 1e-8);
//! ```

use crate::units::{define_unit, Quantity};

use crate::units::{Dimension, Unit};
pub enum Length {}
impl Dimension for Length {}
pub trait LengthUnit: Unit<Dim = Length> {}
impl<T: Unit<Dim = Length>> LengthUnit for T {}


define_unit!("m", Meter, Length, 1.0);
pub type Meters = Quantity<Meter>;

define_unit!("Km", Kilometer, Length, 1_000.0);
pub type Km = Kilometer;
pub type Kilometers = Quantity<Km>;
pub const KM: Kilometers = Kilometers::new(1.0);

define_unit!("Au", AstronomicalUnit, Length, 149_597_870.7);
pub type Au = AstronomicalUnit;
pub type AstronomicalUnits = Quantity<Au>;
pub const AU: AstronomicalUnits = AstronomicalUnits::new(1.0);

define_unit!("Ly", LightYear, Length, 9_460_730_472_580.8);
pub type Ly = LightYear;
pub type LightYears = Quantity<Ly>;
pub const LY: LightYears = LightYears::new(1.0);

define_unit!("SR", SolarRadius, Length, 695_700_000.0);
pub type SolarRadiuses = Quantity<SolarRadius>;
pub const SR: SolarRadiuses = SolarRadiuses::new(1.0);


/// AstronomicalUnit -> LightYear.
impl From<AstronomicalUnits> for LightYears { fn from(au: AstronomicalUnits) -> Self { au.to::<LightYear>() } }

/// LightYear -> AstronomicalUnits.
impl From<LightYears> for AstronomicalUnits { fn from(ly: LightYears) -> Self { ly.to::<AstronomicalUnit>() } }

// TODO: Remove me
impl Unit for f64 {
    const RATIO: f64 = 1.0;
    type Dim = Length;
    const SYMBOL: &'static str = "";
}
