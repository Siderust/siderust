// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Asteroid module
//!
//! ## Scientific scope
//!
//! This module models **asteroids / minor planets** — small Solar System
//! bodies orbiting the Sun, primarily in the main belt between Mars and
//! Jupiter or in near-Earth space.  Each body carries a formal designation,
//! taxonomic composition class, Keplerian orbital elements, and an optional
//! Bond albedo.
//!
//! | Constant | Designation | Class | Potentially Hazardous? | Notes |
//! |----------|-------------|-------|-----------------------|-------|
//! | [`CERES_AST`] | (1) Ceres | Main‑belt / Dwarf‑planet | No | Largest asteroid, also dwarf planet |
//! | [`BENNU`] | (101955) Bennu | Near‑Earth (Apollo) | **Yes** | OSIRIS‑REx sample‑return target |
//! | [`APOPHIS`] | (99942) Apophis | Near‑Earth (Aten) | **Yes** | 2029 Earth fly‑by ~0.1 LD |
//!
//! ## Technical scope
//!
//! - [`Asteroid`] — runtime struct with lifetimed string fields and an optional
//!   Bond albedo ([`Albedos`]).
//! - [`AsteroidClass`] — taxonomic enum (`MainBelt`, `NearEarth`, …).
//! - [`AsteroidBuilder`] — builder for runtime construction.
//! - [`Asteroid::with_albedo`] — `const`-safe albedo builder method.
//!
//! ## References
//!
//! - JPL Small-Body Database. <https://ssd.jpl.nasa.gov/tools/sbdb_query.html>
//! - NASA Planetary Fact Sheet. <https://nssdc.gsfc.nasa.gov/planetary/factsheet/>

use crate::astro::orbit::KeplerianOrbit;
use crate::qtty::{Albedos, AstronomicalUnits, Degrees};
use crate::time::JulianDate;

/// Taxonomic class of a small Solar‑System body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsteroidClass {
    MainBelt,
    NearEarth,
    Trojan,
    Centaur,
    TransNeptunian,
    DwarfPlanet,
}

/// Data structure representing an **asteroid / minor planet**.
#[derive(Clone, Debug)]
pub struct Asteroid<'a> {
    pub name: &'a str,
    pub designation: &'a str,
    pub composition: &'a str,
    pub class: AsteroidClass,
    pub orbit: KeplerianOrbit,
    /// Bond albedo (dimensionless, ∈ [0, 1]).  `None` when not catalogued.
    pub albedo: Option<Albedos>,
}

impl<'a> Asteroid<'a> {
    /// Const constructor used by preset constants.
    pub const fn new_const(
        name: &'a str,
        designation: &'a str,
        composition: &'a str,
        class: AsteroidClass,
        orbit: KeplerianOrbit,
    ) -> Self {
        Self {
            name,
            designation,
            composition,
            class,
            orbit,
            albedo: None,
        }
    }

    /// Set the Bond albedo, returning `self` (const-safe builder step).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::asteroid::CERES_AST;
    /// let a = CERES_AST.albedo.unwrap();
    /// assert!((a.value() - 0.09).abs() < 1e-9);
    /// ```
    pub const fn with_albedo(mut self, albedo: Albedos) -> Self {
        self.albedo = Some(albedo);
        self
    }

    /// Fluent builder.
    pub fn builder() -> AsteroidBuilder<'a> {
        AsteroidBuilder::default()
    }
}

// -------------------------------------------------------------------------------------------------
//  Builder
// -------------------------------------------------------------------------------------------------

#[derive(Default, Clone, Debug)]
pub struct AsteroidBuilder<'a> {
    name: Option<&'a str>,
    designation: Option<&'a str>,
    composition: Option<&'a str>,
    class: Option<AsteroidClass>,
    orbit: Option<KeplerianOrbit>,
    albedo: Option<Albedos>,
}

impl<'a> AsteroidBuilder<'a> {
    pub fn name(mut self, name: &'a str) -> Self {
        self.name = Some(name);
        self
    }
    pub fn designation(mut self, des: &'a str) -> Self {
        self.designation = Some(des);
        self
    }
    pub fn composition(mut self, comp: &'a str) -> Self {
        self.composition = Some(comp);
        self
    }
    pub fn class(mut self, class: AsteroidClass) -> Self {
        self.class = Some(class);
        self
    }
    pub fn orbit(mut self, orbit: KeplerianOrbit) -> Self {
        self.orbit = Some(orbit);
        self
    }
    pub fn albedo(mut self, albedo: Albedos) -> Self {
        self.albedo = Some(albedo);
        self
    }

    pub fn build(self) -> Asteroid<'a> {
        Asteroid {
            name: self.name.expect("missing name"),
            designation: self.designation.expect("missing designation"),
            composition: self.composition.unwrap_or("Unknown"),
            class: self.class.unwrap_or(AsteroidClass::MainBelt),
            orbit: self.orbit.expect("missing orbit"),
            albedo: self.albedo,
        }
    }
}

// -------------------------------------------------------------------------------------------------
//  Preset asteroids (orbital elements from JPL SBDB, epoch J2000)
// -------------------------------------------------------------------------------------------------

/// **(1) Ceres** – largest object in the asteroid belt (also classified as a dwarf planet).
pub const CERES_AST: Asteroid = Asteroid::new_const(
    "Ceres",
    "(1) Ceres",
    "C‑type (carbonaceous)",
    AsteroidClass::DwarfPlanet,
    KeplerianOrbit::new(
        AstronomicalUnits::new(2.7675),
        0.0758,
        Degrees::new(10.5941),
        Degrees::new(80.3055),
        Degrees::new(73.5977),
        Degrees::new(95.9892),
        JulianDate::J2000,
    ),
).with_albedo(Albedos::new(0.09));

/// **(101955) Bennu** – OSIRIS‑REx target; carbon‑rich Apollo NEO.
pub const BENNU: Asteroid = Asteroid::new_const(
    "Bennu",
    "(101955) Bennu",
    "B‑type (carbonaceous)",
    AsteroidClass::NearEarth,
    KeplerianOrbit::new(
        AstronomicalUnits::new(1.1264),
        0.203_745,
        Degrees::new(6.03494),
        Degrees::new(2.06084),
        Degrees::new(66.2221),
        Degrees::new(101.703),
        JulianDate::J2000,
    ),
).with_albedo(Albedos::new(0.044));

/// **(99942) Apophis** – Aten NEO with close Earth approach in 2029.
pub const APOPHIS: Asteroid = Asteroid::new_const(
    "Apophis",
    "(99942) Apophis",
    "S‑type (siliceous)",
    AsteroidClass::NearEarth,
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.922495),
        0.191_197,
        Degrees::new(3.33146),
        Degrees::new(204.4722),
        Degrees::new(126.4001),
        Degrees::new(204.479),
        JulianDate::J2000,
    ),
).with_albedo(Albedos::new(0.23));

/// Convenience array of preset asteroids.
pub const ASTEROID_PRESETS: &[&Asteroid] = &[&CERES_AST, &BENNU, &APOPHIS];
