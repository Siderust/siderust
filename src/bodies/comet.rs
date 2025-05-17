//! Comet type and related functionality.
//!
//! Represents comets with their name, tail length, and orbital parameters.
//! - `name`: String identifier for the comet.
//! - `tail_length`: Tail length in kilometers (`Kilometers`).
//! - `orbit`: Orbital parameters (see [`Orbit`]), using SI units (e.g., AU, degrees, Julian Day).

use crate::astro::orbit::Orbit;
use crate::units::Kilometers;

/// Represents a **Comet** characterized by its tail length and orbit.
#[derive(Debug, Clone)]
pub struct Comet {
    pub name: String,
    pub tail_length: Kilometers,
    pub orbit: Orbit,
}

impl Comet {
    pub fn new(name: String, tail_length: Kilometers, orbit: Orbit) -> Self {
        Self {
            name,
            tail_length,
            orbit,
        }
    }
}
