//! Star type and related functionality.
//!
//! Represents stars with name, distance, mass, radius, luminosity, and target.
//! - `name`: Name of the star (borrowed or owned).
//! - `distance`: Distance from Earth in light-years (`LightYear`).
//! - `mass`: Stellar mass in solar masses (`SolarMass`).
//! - `radius`: Stellar radius in solar radii (`SolarRadius`).
//! - `lumminosity`: Stellar luminosity in solar luminosities (`SolarLuminosity`).
//! - `target`: [`Target`] pointing to a Spherical coordinates (see [`SphericalCoord`]), using degrees and Julian Day.

use crate::coordinates::{centers::Geocentric, frames::Equatorial, SphericalCoord};
use crate::targets::Target;
use crate::units::*;

use std::borrow::Cow;

/// Represents a **Star** characterized by its distance, mass, radius, luminosity and position in the sky.
#[derive(Debug, Clone)]
pub struct Star<'a> {
    pub name: Cow<'a, str>,
    pub distance: LightYear,
    pub mass: SolarMass,
    pub radius: SolarRadius,
    pub lumminosity: SolarLuminosity,
    pub target: Target<SphericalCoord<Geocentric, Equatorial>>,
}

impl<'a> Star<'a> {
    /// Compile‐time constructor: only works with `'static` string literals.
    pub const fn new_const(
        name: &'static str,
        distance: LightYear,
        mass: SolarMass,
        radius: SolarRadius,
        lumminosity: SolarLuminosity,
        target: Target<SphericalCoord<Geocentric, Equatorial>>,
    ) -> Star<'static> {
        Star {
            name: Cow::Borrowed(name),
            distance,
            mass,
            radius,
            lumminosity,
            target,
        }
    }

    /// Runtime constructor: accepts any string-like thing.
    pub fn new<N>(
        name: N,
        distance: LightYear,
        mass: SolarMass,
        radius: SolarRadius,
        lumminosity: SolarLuminosity,
        target: Target<SphericalCoord<Geocentric, Equatorial>>,
    ) -> Star<'a>
    where
        N: Into<Cow<'a, str>>,
    {
        Star {
            name: name.into(),
            distance,
            mass,
            radius,
            lumminosity,
            target,
        }
    }
}
