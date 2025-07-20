//! # Astronomical Bodies
//!
//! A **central hub** that exposes type definitions, ephemerides, and physical
//! parameters for the natural objects that populate our Universe.  Whether you
//! need the sidereal period of *Mars*, the absolute magnitude of *Halley’s Comet*,
//! or the radius of *Betelgeuse*, this module places those numbers at your
//! fingertips while preserving **strong type‑safety** through siderust unit system.
//!
//! The module is intentionally lightweight: all frequently accessed constants
//! (e.g. `SUN`, `MOON`, `JUPITER`) are available at **compile‑time** so they can
//! be used inside `const fn` and embedded targets. More detailed or mutable
//! datasets—such as time‑dependent orbital elements—live in their dedicated
//! sub‑modules but are re‑exported here for convenience.
//!
//! ## Sub‑module overview
//! | Path | Purpose | Highlights |
//! |------|---------|------------|
//! | [`solar_system`] | Canonical Sun, planets, major moons, dwarf planets, and Lagrange points | `SOLAR_SYSTEM` aggregate constant |
//! | [`planets`] | Generic [`Planet`] type plus planet builder utilities | Unit‑safe Keplerian elements |
//! | [`satelite`] | Natural moons and artificial satellites (TLE support) | `MOON` constant, SGP4 propagator |
//! | [`comet`] | Barycentric and heliocentric comet orbits | Halley, Encke, Hale‑Bopp presets |
//! | [`asteroid`] | Minor planets and near‑Earth objects | Ceres, Bennu, Apophis |
//! | [`stars`] | Generic [`Star`] type and stellar parameter helpers |  mass, luminosity ... |
//! | [`catalog`] | Hand‑picked bright‑star catalog | [`catalog::VEGA`], [`catalog::SIRIUS`] |
//!
//! ## Quick start
//! ```rust
//! use siderust::bodies::{Star, Planet, EARTH};
//!
//! // Look up Earth in the pre‑baked Solar System constant
//! let earth: &Planet = &EARTH;
//! println!("Earth {:?}", earth);
//!
//! // Grab a famous bright star from the catalog
//! use siderust::bodies::catalog::ALTAIR;
//! println!("{} is {} away", ALTAIR.name, ALTAIR.distance);
//! ```
//!
//! ## Design philosophy
//! * **Type safety first** — every physical quantity carries its unit.
//! * **Zero runtime cost** for built‑in constants.
//! * **Ergonomic**: most commonly used types (`Planet`, `Star`, `Comet`, …) are
//!   re‑exported at the crate root to minimise path noise.
//! * **Extensibility**: each sub‑module can be expanded independently as new
//!   data become available.
//!
//! ## Re‑exports
//! The following items are publicly re‑exported so you can write, for example,
//! `use astro::bodies::Star;` instead of `use astro::bodies::stars::Star;`:
//!
//! * [`Comet`]
//! * [`Asteroid`]
//! * [`Satellite`]
//! * [`Planet`]
//! * [`Star`]
//! * Everything inside [`solar_system`]
//!
//! ---

pub mod asteroid;
pub mod satelite;
pub mod stars;
pub mod comet;
pub mod planets;
pub mod catalog;
pub mod solar_system;

pub use comet::Comet;
pub use asteroid::Asteroid;
pub use satelite::Satellite;
pub use planets::Planet;
pub use stars::Star;

pub use catalog::*;
pub use solar_system::*;
