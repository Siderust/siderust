//! # Mass Units Module
//!
//! This module provides types and utilities for handling mass-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! mass units and conversions between them.
//!
//! ## Features
//! - **Kilograms (kg)**: The SI base unit of mass, with arithmetic operations.
//! - **Solar Masses (M☉)**: The mass of the Sun, commonly used in astronomy.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{Kilograms, SolarMasses};
//!
//! let m_kg = Kilograms::new(5.0);
//! assert_eq!(m_kg.value(), 5.0);
//!
//! let m_solar = SolarMasses::new(2.0);
//! assert_eq!(m_solar.value(), 2.0);
//! ```

use crate::units::{define_unit, Quantity};

use crate::units::{Dimension, Unit};
pub enum Mass {}
impl Dimension for Mass {}
pub trait MassUnit: Unit<Dim = Mass> {}
impl<T: Unit<Dim = Mass>> MassUnit for T {}

define_unit!("g", Gram, Mass, 1.0);
pub type Grams = Quantity<Gram>;

define_unit!("Kg", Kilogram, Mass, 1_000.0);
pub type Kg = Kilogram;
pub type Kilograms = Quantity<Kg>;
pub const KG: Kilograms = Kilograms::new(1.0);

define_unit!("M☉", SolarMass, Mass, 1.98847e33);
pub type SolarMasses = Quantity<SolarMass>;