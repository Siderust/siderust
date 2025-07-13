//! Asteroid type and related functionality.
//!
//! Represents asteroids with their name, composition, and orbital parameters.
//! - `name`: String identifier for the asteroid.
//! - `composition`: String describing the primary material makeup.
//! - `orbit`: Orbital parameters (see [`Orbit`]), using SI units (e.g., AstronomicalUnits, degrees, Julian Day).

use crate::astro::orbit::Orbit;

/// Represents an **asteroid** characterized by its composition and orbit.
#[derive(Clone)]
pub struct Asteroid {
    pub name: String,
    pub composition: String,
    pub orbit: Orbit,
}

impl Asteroid {
    pub fn new(name: String, composition: String, orbit: Orbit) -> Self {
        Self {
            name,
            composition,
            orbit,
        }
    }
}