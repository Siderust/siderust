//! Satellite type and related functionality.
//!
//! Represents natural or artificial satellites with name, mass, radius, and orbital parameters.
//! - `name`: Name of the satellite (borrowed or owned).
//! - `mass`: Mass in kilograms (`Kilograms`).
//! - `radius`: Mean radius in kilometers (`Kilometers`).
//! - `orbit`: Orbital parameters (see [`Orbit`]).

use crate::astro::orbit::Orbit;
use crate::units::{Kilograms, Kilometers};

use std::borrow::Cow;

#[derive(Clone, Debug)]
/// Represents a **Satelite** characterized by its mass, radius and orbit.
pub struct Satelite<'a> {
    pub name: Cow<'a, str>,
    pub mass: Kilograms,
    pub radius: Kilometers,
    pub orbit: Orbit,
}

impl<'a> Satelite<'a> {
    /// Compile‐time constructor: only works with `'static` string literals.
    pub const fn new_const(
        name: &'static str,
        mass: Kilograms,
        radius: Kilometers,
        orbit: Orbit,
    ) -> Satelite<'static> {
        Satelite {
            name: Cow::Borrowed(name),
            mass,
            radius,
            orbit,
        }
    }

    /// Runtime constructor: accepts any string-like thing.
    pub fn new<N>(name: N, mass: Kilograms, radius: Kilometers, orbit: Orbit) -> Satelite<'a>
    where
        N: Into<Cow<'a, str>>,
    {
        Satelite {
            name: name.into(),
            mass,
            radius,
            orbit,
        }
    }
}
