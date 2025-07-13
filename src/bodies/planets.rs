//! Planet type and related functionality.
//!
//! Represents planets with mass, radius, and orbital parameters.
//! - `mass`: Planetary mass in kilograms (`Kilograms`).
//! - `radius`: Mean planetary radius in kilometers (`Kilometers`).
//! - `orbit`: Orbital parameters (see [`Orbit`]), using SI units (e.g., AstronomicalUnits, degrees, Julian Day).

use crate::astro::orbit::Orbit;
use crate::units::{Kilometers, Kilograms};

/// Represents a **Planet** characterized by its mass, radius and orbit.
#[derive(Clone)]
pub struct Planet {
    pub mass: Kilograms,
    pub radius: Kilometers,
    pub orbit: Orbit,
}
