//! Data and helper types for natural and artificial bodies.
//!
//! The `bodies` module exposes ephemerides and physical parameters for a range
//! of solar‑system objects and stars.  Frequently used constants such as
//! `SUN`, `MOON` or `JUPITER` are available as `const` values for use in
//! compile‑time contexts, while more extensive datasets live in dedicated
//! submodules and are re‑exported here for convenience.  All quantities retain
//! their units through the [`units`](../units/index.html) system to avoid
//! accidental mixing.
//!
//! ## Submodules
//! - [`solar_system`]: canonical data for planets, major moons and Lagrange points.
//! - [`planets`]: [`Planet`] type and builders from Keplerian elements.
//! - [`satelite`]: natural satellites and artificial spacecraft with optional SGP4 propagation.
//! - [`comet`]: preset and user‑defined comet orbits.
//! - [`asteroid`]: minor‑planet definitions.
//! - [`stars`]: [`Star`] type and stellar parameters.
//! - [`catalog`]: curated catalog of bright stars.
//!
//! ## Example
//! ```rust
//! use siderust::bodies::{catalog::ALTAIR, EARTH};
//!
//! // Access Earth from the solar‑system constants
//! println!("Earth radius: {}", EARTH.radius);
//!
//! // Retrieve a bright star from the catalog
//! println!("{} is {} away", ALTAIR.name, ALTAIR.distance);
//! ```
//!
//! ## Re‑exports
//! [`Comet`], [`Asteroid`], [`Satellite`], [`Planet`], [`Star`] and all items
//! from [`solar_system`].

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
