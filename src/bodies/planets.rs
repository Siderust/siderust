// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Planets
//!
//! ## Scientific scope
//!
//! Planets in our Solar System range from the rocky inner worlds (Mercury,
//! Venus, Earth, Mars) to the gas and ice giants of the outer system.
//! This module models planets as bodies defined by mass, mean radius, and
//! heliocentric Keplerian orbital elements at J2000.0. The Bond albedo —
//! the fraction of total incident solar irradiance reflected back to space —
//! is an optional typed field that supports photometric modelling.
//!
//! Keplerian orbital periods are derived from Gauss's form of Kepler's
//! third law:
//!
//! ```text
//! P = (2π / k) sqrt(a³),   k ≈ 0.017 202 098 95 AU^{3/2} d^{-1}
//! ```
//!
//! where `a` is the semi-major axis in AU and `P` is returned in seconds.
//!
//! ## Technical scope
//!
//! - [`Planet`] — data structure for mass ([`Kilograms`]), radius
//!   ([`Kilometers`]), orbit ([`KeplerianOrbit`]), and optional Bond
//!   albedo ([`Albedos`]).
//! - [`PlanetBuilder`] / [`PlanetBuilderError`] — fluent builder with
//!   explicit error variants.
//! - [`OrbitExt`] — extension trait on [`KeplerianOrbit`] providing
//!   [`OrbitExt::period`] (typed [`Seconds`]).
//!
//! Planet constants (Mercury … Neptune, dwarf planets, moons) live in
//! [`crate::bodies::solar_system`]; this module only provides the generic
//! type and builder.
//!
//! ## References
//!
//! - Gauss, C. F. (1809). *Theoria Motus Corporum Coelestium*. Hamburg.
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. Willmann-Bell.
//!   Chapter 33.
//! - IAU Working Group on Cartographic Coordinates and Rotational Elements
//!   (2015). *Celestial Mechanics and Dynamical Astronomy* 130, 22.
//!   doi:10.1007/s10569-017-9805-5
//!
//! ## Quick start
//! ```rust
//! use siderust::bodies::planets::{Planet, PlanetBuilder};
//! use siderust::qtty::*;
//! use siderust::astro::orbit::KeplerianOrbit;
//! use siderust::time::JulianDate;
//!
//! let kepler_452b = Planet::builder()
//!     .mass(Kilograms::new(5.972e24 * 5.0))
//!     .radius(Kilometers::new(6371.0 * 1.6))
//!     .orbit(KeplerianOrbit::new(
//!         AstronomicalUnits::new(1.046),
//!         0.02,
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         Degrees::new(0.0),
//!         siderust::J2000,
//!     ))
//!     .build();
//!
//! println!("Custom planet mass: {}", kepler_452b.mass);
//! ```

use crate::astro::orbit::KeplerianOrbit;
use crate::qtty::*;

/// Represents a **Planet** characterised by its mass, mean radius, Bond albedo
/// and Keplerian heliocentric orbit.
#[derive(Clone, Debug)]
pub struct Planet {
    pub mass: Kilograms,
    pub radius: Kilometers,
    pub orbit: KeplerianOrbit,
    /// Bond albedo (dimensionless, ∈ [0, 1]).  `None` when not catalogued.
    pub albedo: Option<Albedos>,
}

impl Planet {
    /// Compile-time constructor (`albedo` defaults to `None`).
    ///
    /// Use [`Planet::with_albedo`] in a `const` context to attach an albedo.
    pub const fn new_const(mass: Kilograms, radius: Kilometers, orbit: KeplerianOrbit) -> Self {
        Self {
            mass,
            radius,
            orbit,
            albedo: None,
        }
    }

    /// Attach a typed Bond albedo (`const`-safe builder).
    ///
    /// # Arguments
    ///
    /// - `albedo` — Bond albedo as [`Albedos`] (`Quantity<Albedo>`), ∈ [0, 1].
    pub const fn with_albedo(mut self, albedo: Albedos) -> Self {
        self.albedo = Some(albedo);
        self
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
    orbit: Option<KeplerianOrbit>,
    albedo: Option<Albedos>,
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
    pub fn orbit(mut self, orbit: KeplerianOrbit) -> Self {
        self.orbit = Some(orbit);
        self
    }

    /// Set Bond albedo (dimensionless, ∈ [0, 1]).
    pub fn albedo(mut self, albedo: Albedos) -> Self {
        self.albedo = Some(albedo);
        self
    }

    /// Attempt to build a [`Planet`], returning an error if any field is missing.
    pub fn try_build(self) -> Result<Planet, PlanetBuilderError> {
        Ok(Planet {
            mass: self.mass.ok_or(PlanetBuilderError::MissingMass)?,
            radius: self.radius.ok_or(PlanetBuilderError::MissingRadius)?,
            orbit: self.orbit.ok_or(PlanetBuilderError::MissingOrbit)?,
            albedo: self.albedo,
        })
    }

    /// Build a [`Planet`], panicking on missing fields.  Suitable for tests
    /// and quick demos where construction completeness is guaranteed.
    pub fn build(self) -> Planet {
        self.try_build().expect("incomplete PlanetBuilder")
    }

    /// Const-context helper when all three mandatory fields are already present.
    pub const fn build_unchecked(self) -> Planet {
        match (self.mass, self.radius, self.orbit) {
            (Some(mass), Some(radius), Some(orbit)) => Planet {
                mass,
                radius,
                orbit,
                albedo: None,
            },
            _ => panic!("PlanetBuilder::build_unchecked called with missing fields"),
        }
    }
}

// -------------------------------------------------------------------------------------------------
//  Extension trait for derived orbital quantities
// -------------------------------------------------------------------------------------------------

/// Additional derived methods for [`KeplerianOrbit`].
pub trait OrbitExt {
    /// Sidereal orbital period using Kepler's third law, returned in
    /// typed [`Seconds`].
    ///
    /// # Returns
    ///
    /// Period in [`Seconds`].
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::{EARTH, planets::OrbitExt};
    ///
    /// let p = EARTH.orbit.period();
    /// // Earth period ≈ 365.25 days = 31 557 600 s
    /// assert!((p.value() - 3.156e7).abs() < 1e5);
    /// ```
    fn period(&self) -> Seconds;
}

impl OrbitExt for KeplerianOrbit {
    fn period(&self) -> Seconds {
        let a_au = self
            .shape()
            .semi_major_axis()
            .to::<AstronomicalUnit>()
            .value();
        Seconds::new(crate::astro::units::heliocentric_period_days(a_au) * 86_400.0)
    }
}

// -------------------------------------------------------------------------------------------------
//  Unit tests
// -------------------------------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{AstronomicalUnits, Degrees, Kilograms, Kilometers};

    #[test]
    fn builder_roundtrip() {
        let p = Planet::builder()
            .mass(Kilograms::new(1.0))
            .radius(Kilometers::new(1.0))
            .orbit(KeplerianOrbit::new(
                AstronomicalUnits::new(1.0),
                0.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                crate::J2000,
            ))
            .build();
        assert_eq!(p.mass.value(), 1.0);
    }

    #[test]
    fn with_albedo_roundtrip() {
        let p = Planet::builder()
            .mass(Kilograms::new(1.0))
            .radius(Kilometers::new(1.0))
            .orbit(KeplerianOrbit::new(
                AstronomicalUnits::new(1.0),
                0.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                crate::J2000,
            ))
            .albedo(Albedos::new(0.30))
            .build();
        assert!((p.albedo.unwrap().value() - 0.30).abs() < 1e-10);
    }
}
