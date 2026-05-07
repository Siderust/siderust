// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Satellites
//!
//! ## Scientific scope
//!
//! Natural satellites (moons) orbit planets or dwarf planets under
//! gravitational attraction. Their orbits are often non-Keplerian due to
//! oblateness of the primary and mutual perturbations, but for many
//! practical purposes a mean Keplerian element set gives adequate
//! first-order positions. The Bond albedo characterises the overall
//! reflectivity of the body's surface integrated over all wavelengths and
//! phase angles.
//!
//! This module does **not** model artificial satellites; TLE-based SGP4
//! propagation is outside scope here.
//!
//! ## Technical scope
//!
//! - [`Satellite`] — data structure for name, mass ([`Kilograms`]),
//!   mean radius ([`Kilometers`]), Keplerian orbit ([`KeplerianOrbit`]),
//!   and optional Bond albedo ([`Albedos`]).
//! - [`Satellite::new_const`] — `const`-safe constructor.
//! - [`Satellite::with_albedo`] — `const`-safe albedo builder.
//! - [`Satellite::new`] — runtime constructor accepting any string-like name.
//!
//! Pre-built satellite constants (Moon, Io, Europa, Ganymede, Callisto,
//! Titan, Triton) are defined in [`crate::bodies::solar_system`].
//!
//! ## References
//!
//! - Seidelmann, P. K. (Ed.) (1992). *Explanatory Supplement to the
//!   Astronomical Almanac*. University Science Books. Chapter 15.
//! - IAU Working Group on Cartographic Coordinates and Rotational Elements
//!   (2015). *Celestial Mechanics and Dynamical Astronomy* 130, 22.
//!   doi:10.1007/s10569-017-9805-5

use crate::astro::orbit::KeplerianOrbit;
use crate::qtty::{Albedos, Kilograms, Kilometers};

use std::borrow::Cow;

/// Represents a **Satellite** characterized by its mass, radius, optional
/// Bond albedo, and Keplerian orbit.
#[derive(Clone, Debug)]
pub struct Satellite<'a> {
    pub name: Cow<'a, str>,
    pub mass: Kilograms,
    pub radius: Kilometers,
    pub orbit: KeplerianOrbit,
    /// Bond albedo (dimensionless, ∈ [0, 1]).  `None` when not catalogued.
    pub albedo: Option<Albedos>,
}

impl<'a> Satellite<'a> {
    /// Compile-time constructor (only works with `'static` string literals).
    ///
    /// `albedo` defaults to `None`; use [`Satellite::with_albedo`] to
    /// attach a known Bond albedo value.
    pub const fn new_const(
        name: &'static str,
        mass: Kilograms,
        radius: Kilometers,
        orbit: KeplerianOrbit,
    ) -> Satellite<'static> {
        Satellite {
            name: Cow::Borrowed(name),
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

    /// Runtime constructor: accepts any string-like thing.
    pub fn new<N>(
        name: N,
        mass: Kilograms,
        radius: Kilometers,
        orbit: KeplerianOrbit,
    ) -> Satellite<'a>
    where
        N: Into<Cow<'a, str>>,
    {
        Satellite {
            name: name.into(),
            mass,
            radius,
            orbit,
            albedo: None,
        }
    }
}
