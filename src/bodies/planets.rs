//! # Planet module
//!
//! Provides a **generic [`Planet`] type** together with ergonomic *builder*
//! utilities and helper traits for working with unit‑safe Keplerian elements
//! ([`Orbit`]).  These additions let you construct planetary bodies in a
//! progressively‑filled, compile‑time‑safe manner while keeping the original
//! lightweight `struct` used by the Solar‑System constants fully compatible.
//!
//! ## Quick start
//! ```rust
//! use siderust::bodies::planets::{Planet, PlanetBuilder};
//! use siderust::units::{Kilograms, Kilometers, AstronomicalUnits, Degrees};
//! use siderust::astro::orbit::Orbit;
//! use siderust::astro::JulianDate;
//!
//! // Build a custom exoplanet in two steps:
//! let kepler_452b = Planet::builder()
//!     .mass(Kilograms::new(5.972e24 * 5.0))
//!     .radius(Kilometers::new(6371.0 * 1.6))
//!     .orbit(Orbit::new(
//!         AstronomicalUnits::new(1.046),
//!         0.02,
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         JulianDate::J2000,
//!     ))
//!     .build();
//!
//! println!("Custom planet mass: {}", kepler_452b.mass);
//! ```
//!
//! ## Builder pattern
//! The [`PlanetBuilder`] collects optional fields (`mass`, `radius`, `orbit`) and
//! validates their presence when `.build()` is invoked.  Missing fields trigger
//! a [`PlanetBuilderError`]. For convenience, `.build_unchecked()` exists when
//! the caller can guarantee completeness at compile time (e.g. const‑context
//! constructions).
//!
//! ## Helper traits
//! A small extension trait [`OrbitExt`] is provided to compute derived orbital
//! quantities—currently the sidereal period via Kepler’s third law. The trait is
//! implemented for [`Orbit`], so any existing orbit constant gets the new
//! methods *for free*.
//! ---

use crate::astro::orbit::Orbit;
use crate::units::*;

const GAUSSIAN_GRAVITATIONAL_CONSTANT: RadiansPerDay = RadiansPerDay::new(0.017_202_098_95);

/// Represents a **Planet** characterised by its mass, mean radius, and orbit.
#[derive(Clone, Debug)]
pub struct Planet {
    pub mass: Kilograms,
    pub radius: Kilometers,
    pub orbit: Orbit,
}

impl Planet {
    /// Compile‑time constructor.
    pub const fn new_const(mass: Kilograms, radius: Kilometers, orbit: Orbit) -> Self {
        Self {
            mass,
            radius,
            orbit,
        }
    }

    /// Start building a planet with the fluent [`PlanetBuilder`] API.
    pub fn builder() -> PlanetBuilder {
        PlanetBuilder::default()
    }
}

// -------------------------------------------------------------------------------------------------
//  Builder implementation
// -------------------------------------------------------------------------------------------------

/// Error returned when mandatory fields are missing in [`PlanetBuilder`].
#[derive(Debug, Clone)]
pub enum PlanetBuilderError {
    MissingMass,
    MissingRadius,
    MissingOrbit,
}

/// Fluent builder for [`Planet`].  All setters accept any type that can be
/// converted into the required unit wrapper via `Into`.
#[derive(Debug, Default, Clone)]
pub struct PlanetBuilder {
    mass: Option<Kilograms>,
    radius: Option<Kilometers>,
    orbit: Option<Orbit>,
}

impl PlanetBuilder {
    /// Set planetary mass.
    pub fn mass(mut self, mass: impl Into<Kilograms>) -> Self {
        self.mass = Some(mass.into());
        self
    }

    /// Set mean planetary radius.
    pub fn radius(mut self, radius: impl Into<Kilometers>) -> Self {
        self.radius = Some(radius.into());
        self
    }

    /// Set the heliocentric orbit.
    pub fn orbit(mut self, orbit: Orbit) -> Self {
        self.orbit = Some(orbit);
        self
    }

    /// Attempt to build a [`Planet`], returning an error if any field is missing.
    pub fn try_build(self) -> Result<Planet, PlanetBuilderError> {
        Ok(Planet {
            mass: self.mass.ok_or(PlanetBuilderError::MissingMass)?,
            radius: self.radius.ok_or(PlanetBuilderError::MissingRadius)?,
            orbit: self.orbit.ok_or(PlanetBuilderError::MissingOrbit)?,
        })
    }

    /// Build a [`Planet`], panicking on missing fields.  Suitable for tests and
    /// quick demos where construction completeness is guaranteed.
    pub fn build(self) -> Planet {
        self.try_build().expect("incomplete PlanetBuilder")
    }

    /// Const‑context helper when all three fields are already present.
    pub const fn build_unchecked(self) -> Planet {
        match (self.mass, self.radius, self.orbit) {
            (Some(mass), Some(radius), Some(orbit)) => Planet {
                mass,
                radius,
                orbit,
            },
            _ => panic!("PlanetBuilder::build_unchecked called with missing fields"),
        }
    }
}

// -------------------------------------------------------------------------------------------------
//  Extension trait for derived orbital quantities
// -------------------------------------------------------------------------------------------------

/// Additional derived methods for [`Orbit`].
pub trait OrbitExt {
    /// Sidereal orbital period using Kepler’s third law, returned in **seconds**.
    fn period(&self) -> Seconds;
}

impl OrbitExt for Orbit {
    fn period(&self) -> Seconds {
        // Kepler's 3rd: P = 2π * sqrt(a^3 / μ)
        // Using Gaussian gravitational constant k ≈ 0.01720209895 AU^{3/2} d^{-1}
        // ⇒ period (days) = 2π / k * sqrt(a^3)
        // We compute in seconds directly.

        use std::f64::consts::PI;
        let a_au = self.semi_major_axis.to::<AstronomicalUnit>().value();
        let k = GAUSSIAN_GRAVITATIONAL_CONSTANT.to::<RadianPerDay>().value();

        let t_days = (2.0 * PI / k) * a_au.powf(1.5);
        Seconds::new(t_days * 86_400.0)
    }
}

// -------------------------------------------------------------------------------------------------
//  Unit tests
// -------------------------------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::JulianDate;
    use crate::units::{AstronomicalUnits, Degrees, Kilograms, Kilometers};

    #[test]
    fn builder_roundtrip() {
        let p = Planet::builder()
            .mass(Kilograms::new(1.0))
            .radius(Kilometers::new(1.0))
            .orbit(Orbit::new(
                AstronomicalUnits::new(1.0),
                0.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            ))
            .build();
        assert_eq!(p.mass.value(), 1.0);
    }
}
