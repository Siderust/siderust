//! Astronomical bodies module.
//!
//! This module provides types and data for representing various astronomical bodies,
//! including stars, planets, satellites, comets, asteroids, and the solar system.
//! It exposes both data structures and constants for well-known objects such as the Sun, Moon,
//! planets, and notable stars. The module is organized into submodules for each category of body.

mod comet;
mod asteroid;
mod satelite;
mod planets;
mod stars;
pub mod catalog;
pub mod solar_system;

pub use comet::Comet;
pub use asteroid::Asteroid;
pub use satelite::Satelite;
pub use planets::Planet;
pub use stars::Star;

pub use solar_system::*;
