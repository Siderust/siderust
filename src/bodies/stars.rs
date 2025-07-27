//! Star type and related functionality.
//!
//! Represents stars with name, distance, mass, radius, luminosity, and target.
//! - `name`: Name of the star (borrowed or owned).
//! - `distance`: LengthUnit from Earth in light-years (`LightYear`).
//! - `mass`: Stellar mass in solar masses (`SolarMasses`).
//! - `radius`: Stellar radius in solar radii (`SolarRadiuses`).
//! - `luminosity`: Stellar luminosity in solar luminosities (`SolarLuminosity`).
//! - `target`: [`Target`] pointing to a Spherical coordinates (see [`Position`]), using degrees and Julian Day.

use crate::coordinates::{centers::Geocentric, frames::Equatorial, spherical::Position};
use crate::targets::Target;
use crate::units::*;

use std::borrow::Cow;

/// Represents a **Star** characterized by its distance, mass, radius, luminosity and position in the sky.
#[derive(Clone, Debug)]
pub struct Star<'a> {
    pub name: Cow<'a, str>,
    pub distance: LightYears,
    pub mass: SolarMasses,
    pub radius: SolarRadiuses,
    pub luminosity: SolarLuminosities,
    pub target: Target<Position<Geocentric, Equatorial, LightYear>>,
}

impl<'a> Star<'a> {
    /// Compile‐time constructor: only works with `'static` string literals.
    pub const fn new_const(
        name: &'static str,
        distance: LightYears,
        mass: SolarMasses,
        radius: SolarRadiuses,
        luminosity: SolarLuminosities,
        target: Target<Position<Geocentric, Equatorial, LightYear>>,
    ) -> Star<'static> {
        Star {
            name: Cow::Borrowed(name),
            distance,
            mass,
            radius,
            luminosity,
            target,
        }
    }

    /// Runtime constructor: accepts any string-like thing.
    pub fn new<N>(
        name: N,
        distance: LightYears,
        mass: SolarMasses,
        radius: SolarRadiuses,
        luminosity: SolarLuminosities,
        target: Target<Position<Geocentric, Equatorial, LightYear>>,
    ) -> Star<'a>
    where
        N: Into<Cow<'a, str>>,
    {
        Star {
            name: name.into(),
            distance,
            mass,
            radius,
            luminosity,
            target,
        }
    }
}
