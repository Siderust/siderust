//! # Length Units Module
//!
//! This module provides types and utilities for handling length-related calculations
//! in astronomical and scientific contexts. It includes representations for various
//! length units and conversions between them.
//!
//! ## Features
//! - **Astronomical Unit (AU)**: The mean distance between the Earth and the Sun.
//! - **Light Year (LY)**: The distance light travels in one Julian year in vacuum.
//! - Conversion between AU and LY.
//!
//! ## Example Usage
//! ```rust
//! use siderust::units::{AstronomicalUnit, LightYear};
//!
//! let au = AstronomicalUnit::new(1.0);
//! let ly = LightYear::from(au);
//! assert!((ly.value() - 1.582e-5).abs() < 1e-8);
//! ```

mod astronomical_unit;
mod light_year;
mod kilometer;
mod solar_radius;

pub use astronomical_unit::*;
pub use light_year::*;
pub use kilometer::*;
pub use solar_radius::*;
