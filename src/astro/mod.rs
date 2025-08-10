//! Astronomical time scales and positional corrections.
//!
//! The `astro` module gathers low‑level routines that model phenomena affecting
//! precise astronomical observations.  It covers time representations, Earth
//! orientation, relativistic effects and basic orbital tools.  Each component is
//! allocation‑free and designed for cross‑validation against published models.
//!
//! ## Components
//! - [`julian_date`] and [`modified_julian_date`] time tags commonly used in
//!   dynamical astronomy.
//! - [`dynamical_time`] conversions between universal and terrestrial time.
//! - [`nutation`] and [`precession`] describing Earth orientation.
//! - [`aberration`] annual aberration corrections.
//! - [`sidereal`] computation of Greenwich and local sidereal time.
//! - [`orbit`] helpers for classical orbital elements.
//! - [`proper_motion`] linear motion of stars.

pub mod aberration;
pub mod proper_motion;
pub mod nutation;
pub mod precession;
pub mod sidereal;
pub mod dynamical_time;
pub mod orbit;
pub mod julian_date;
pub mod modified_julian_date;

pub use julian_date::*;
pub use modified_julian_date::*;
