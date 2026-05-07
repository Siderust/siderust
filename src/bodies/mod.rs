// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical Bodies
//!
//! ## Scientific scope
//!
//! This module is the central registry for natural bodies that populate our
//! Solar System and beyond: the Sun, planets, dwarf planets, moons, comets,
//! asteroids, and notable stars. Physical parameters (mass, radius, albedo)
//! and Keplerian orbital elements at J2000.0 are exposed as compile-time
//! constants using strongly-typed unit wrappers from [`crate::qtty`], so
//! dimensional errors are caught at compile time.
//!
//! ## Technical scope
//!
//! This module is a lightweight hub that re-exports the most commonly
//! accessed types and constants.  All frequently used constants (e.g.
//! [`SUN`], [`EARTH`], [`MOON`], [`JUPITER`]) are available at
//! **compile-time** so they can be used inside `const fn` and on
//! embedded targets.  More detailed or mutable datasets (time-dependent
//! orbital elements) live in sub-modules but are re-exported here.
//!
//! | Sub-module | Purpose | Highlights |
//! |------------|---------|------------|
//! | [`solar_system`] | Sun, planets, moons, dwarf planets, Lagrange points | [`SOLAR_SYSTEM`] aggregate |
//! | [`planets`] | Generic [`Planet`] type + builder | Typed Keplerian elements, albedo |
//! | [`satelite`] | Natural moons + artificial satellites | [`MOON`] constant |
//! | [`comet`] | Barycentric and heliocentric comet orbits | [`comet::HALLEY`], [`comet::ENCKE`], [`comet::HALE_BOPP`] |
//! | [`asteroid`] | Minor planets and NEOs | [`asteroid::CERES_AST`], [`asteroid::BENNU`], [`asteroid::APOPHIS`] |
//! | [`stars`] | Generic [`Star`] type and stellar parameter helpers | mass, luminosity |
//! | [`catalog`] | Hand-picked bright-star catalog | [`catalog::VEGA`], [`catalog::SIRIUS`] |
//!
//! ## References
//!
//! - Williams, D. R. (2024). *Planetary Fact Sheet – Metric*. NASA GSFC.
//!   <https://nssdc.gsfc.nasa.gov/planetary/factsheet/>
//! - JPL Small-Body Database (2024). <https://ssd.jpl.nasa.gov/>
//! - IAU (2012). *Resolutions of the 28th General Assembly*. Resolution B2.
//!
//! ## Quick start
//! ```rust
//! use siderust::bodies::{Star, Planet, EARTH};
//!
//! let earth: &Planet = &EARTH;
//! println!("Earth {:?}", earth);
//!
//! use siderust::bodies::catalog::ALTAIR;
//! println!("{} is {} away", ALTAIR.name, ALTAIR.distance);
//! ```

pub mod asteroid;
pub mod catalog;
pub mod comet;
pub mod planets;
pub mod satelite;
pub mod solar_system;
pub mod stars;

pub use asteroid::Asteroid;
pub use comet::Comet;
pub use planets::Planet;
pub use satelite::Satellite;
pub use stars::Star;

pub use catalog::*;
pub use solar_system::*;
