//! Spherical coordinate types.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)

pub mod position;
pub use position::Position;

pub mod direction;
pub use direction::Direction;

mod ecef;
mod ecliptic;
mod equatorial;
mod horizontal;
mod icrs;
