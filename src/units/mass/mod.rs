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
//! use siderust::units::{Kilograms, SolarMass};
//!
//! let m_kg = Kilograms::new(5.0);
//! assert_eq!(m_kg.value(), 5.0);
//!
//! let m_solar = SolarMass::new(2.0);
//! assert_eq!(m_solar.value(), 2.0);
//! ```

mod kilogram;
mod solar_mass;

pub use kilogram::*;
pub use solar_mass::*;
