//! # Asteroid module
//!
//! Represents asteroids with:
//! - `name`: String identifier for the asteroid.
//! - `designation`: system used to formally identify asteroids.
//! - `composition`: String describing the primary material makeup.
//! - `orbit`: Orbital parameters (see [`Orbit`]),.
//!
//! | Constant | Designation | Class | Potentially Hazardous? | Notes |
//! |----------|-------------|-------|-----------------------|-------|
//! | [`CERES_AST`] | (1) Ceres | Main‑belt / Dwarf‑planet | No | Largest asteroid, also dwarf planet |
//! | [`BENNU`] | (101955) Bennu | Near‑Earth ( Apollo ) | **Yes** | OSIRIS‑REx sample‑return target |
//! | [`APOPHIS`] | (99942) Apophis | Near‑Earth ( Aten ) | **Yes** | 2029 Earth fly‑by ~0.1 LD |
//!
//! ## API highlights
//! * `AsteroidClass` enum (`MainBelt`, `NearEarth`, …)
//! * `AsteroidBuilder` for fluent construction of custom objects
//!
//! ---

use crate::astro::orbit::Orbit;
use crate::astro::JulianDate;
use crate::units::{AstronomicalUnits, Degrees};

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
    pub orbit: Orbit,
}

impl<'a> Asteroid<'a> {
    /// Const constructor used by preset constants.
    pub const fn new_const(
        name: &'a str,
        designation: &'a str,
        composition: &'a str,
        class: AsteroidClass,
        orbit: Orbit,
    ) -> Self {
        Self { name, designation, composition, class, orbit }
    }

    /// Fluent builder.
    pub fn builder() -> AsteroidBuilder<'a> { AsteroidBuilder::default() }

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
    orbit: Option<Orbit>,
}

impl<'a> AsteroidBuilder<'a> {
    pub fn name(mut self, name: &'a str) -> Self { self.name = Some(name); self }
    pub fn designation(mut self, des: &'a str) -> Self { self.designation = Some(des); self }
    pub fn composition(mut self, comp: &'a str) -> Self { self.composition = Some(comp); self }
    pub fn class(mut self, class: AsteroidClass) -> Self { self.class = Some(class); self }
    pub fn orbit(mut self, orbit: Orbit) -> Self { self.orbit = Some(orbit); self }

    pub fn build(self) -> Asteroid<'a> {
        Asteroid {
            name: self.name.expect("missing name"),
            designation: self.designation.expect("missing designation"),
            composition: self.composition.unwrap_or("Unknown"),
            class: self.class.unwrap_or(AsteroidClass::MainBelt),
            orbit: self.orbit.expect("missing orbit"),
        }
    }
}

// -------------------------------------------------------------------------------------------------
//  Preset asteroids (orbital elements from JPL SBDB, epoch J2000)
// -------------------------------------------------------------------------------------------------

/// **(1) Ceres** – largest object in the asteroid belt (also classified as a dwarf planet).
pub const CERES_AST: Asteroid = Asteroid::new_const(
    "Ceres",
    "(1) Ceres",
    "C‑type (carbonaceous)",
    AsteroidClass::DwarfPlanet,
    Orbit::new(
        AstronomicalUnits::new(2.7675),
        0.0758,
        Degrees::new(10.5941),
        Degrees::new(80.3055),
        Degrees::new(73.5977),
        Degrees::new(95.9892),
        JulianDate::J2000,
    ),
);

/// **(101955) Bennu** – OSIRIS‑REx target; carbon‑rich Apollo NEO.
pub const BENNU: Asteroid = Asteroid::new_const(
    "Bennu",
    "(101955) Bennu",
    "B‑type (carbonaceous)",
    AsteroidClass::NearEarth,
    Orbit::new(
        AstronomicalUnits::new(1.1264),
        0.203_745,
        Degrees::new(6.03494),
        Degrees::new(2.06084),
        Degrees::new(66.2221),
        Degrees::new(101.703),
        JulianDate::J2000,
    ),
);

/// **(99942) Apophis** – Aten NEO with close Earth approach in 2029.
pub const APOPHIS: Asteroid = Asteroid::new_const(
    "Apophis",
    "(99942) Apophis",
    "S‑type (siliceous)",
    AsteroidClass::NearEarth,
    Orbit::new(
        AstronomicalUnits::new(0.922495),
        0.191_197,
        Degrees::new(3.33146),
        Degrees::new(204.4722),
        Degrees::new(126.4001),
        Degrees::new(204.479),
        JulianDate::J2000,
    ),
);

/// Convenience array of preset asteroids.
pub const ASTEROID_PRESETS: &[&Asteroid] = &[&CERES_AST, &BENNU, &APOPHIS];
