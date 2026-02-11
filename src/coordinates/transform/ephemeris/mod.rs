// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Ephemeris Abstraction Layer
//!
//! This module provides traits for querying planetary positions and velocities
//! from different ephemeris backends (VSOP87, JPL DE, etc.), decoupling the
//! coordinate transformation system from specific ephemeris implementations.
//!
//! ## Design Philosophy
//!
//! The ephemeris layer separates three concerns:
//! - **Geometry**: coordinate translations/rotations (handled by transform module)
//! - **Ephemerides**: time-dependent body state lookup (this module)
//! - **Observation effects**: aberration, refraction, etc. (handled by observation module)
//!
//! ## Usage
//!
//! For most users, the [`Vsop87Ephemeris`] default is sufficient:
//!
//! ```rust
//! use siderust::coordinates::transform::ephemeris::{BodyId, BodyEphemeris, Vsop87Ephemeris};
//! use siderust::astro::JulianDate;
//!
//! let eph = Vsop87Ephemeris;
//! let earth_pos = eph.position_barycentric(BodyId::Earth, JulianDate::J2000);
//! ```
//!
//! ## Implementing Custom Ephemerides
//!
//! To use a different ephemeris source (e.g., JPL DE), implement the
//! [`BodyEphemeris`] trait:
//!
//! ```rust,ignore
//! struct JplDeEphemeris { /* loaded DE data */ }
//!
//! impl BodyEphemeris for JplDeEphemeris {
//!     fn position_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
//!         // Look up position from DE tables
//!     }
//! }
//! ```

mod vsop87;

pub use vsop87::Vsop87Ephemeris;

use crate::astro::JulianDate;

/// Identifies a solar system body for ephemeris queries.
///
/// This enum allows ephemeris providers to handle different bodies
/// through a single interface rather than separate methods per body.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BodyId {
    /// The Sun (center of heliocentric coordinates)
    Sun,
    /// Earth (center of geocentric coordinates)
    Earth,
    /// The Moon (Earth's natural satellite)
    Moon,
    /// Mercury
    Mercury,
    /// Venus
    Venus,
    /// Mars
    Mars,
    /// Jupiter
    Jupiter,
    /// Saturn
    Saturn,
    /// Uranus
    Uranus,
    /// Neptune
    Neptune,
}

impl std::fmt::Display for BodyId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BodyId::Sun => write!(f, "Sun"),
            BodyId::Earth => write!(f, "Earth"),
            BodyId::Moon => write!(f, "Moon"),
            BodyId::Mercury => write!(f, "Mercury"),
            BodyId::Venus => write!(f, "Venus"),
            BodyId::Mars => write!(f, "Mars"),
            BodyId::Jupiter => write!(f, "Jupiter"),
            BodyId::Saturn => write!(f, "Saturn"),
            BodyId::Uranus => write!(f, "Uranus"),
            BodyId::Neptune => write!(f, "Neptune"),
        }
    }
}

/// Trait for querying body positions from an ephemeris source.
///
/// Implementations provide time-dependent positions of solar system bodies
/// in barycentric and/or heliocentric reference frames.
///
/// # Coordinate Convention
///
/// All positions are returned in **ecliptic coordinates** (J2000.0 mean ecliptic)
/// with units of **astronomical units (AU)**.
///
/// # Required Methods
///
/// Implementors must provide at least [`position_barycentric`](Self::position_barycentric).
/// A default implementation of [`position_heliocentric`](Self::position_heliocentric) is
/// provided that computes `helio = bary - sun_bary`.
pub trait BodyEphemeris {
    /// Returns the position of `body` relative to the solar system barycenter.
    ///
    /// # Arguments
    ///
    /// - `body`: The body to query.
    /// - `jd`: The Julian Date (TDB scale preferred, TT acceptable for most uses).
    ///
    /// # Returns
    ///
    /// Position `[x, y, z]` in AU, in ecliptic J2000.0 coordinates.
    fn position_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3];

    /// Returns the position of `body` relative to the Sun (heliocentric).
    ///
    /// # Default Implementation
    ///
    /// Computes `position_barycentric(body) - position_barycentric(Sun)`.
    fn position_heliocentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        let [bx, by, bz] = self.position_barycentric(body, jd);
        let [sx, sy, sz] = self.position_barycentric(BodyId::Sun, jd);
        [bx - sx, by - sy, bz - sz]
    }
}

/// Trait for querying body velocities from an ephemeris source.
///
/// Extends [`BodyEphemeris`] with velocity information, required for
/// computing aberration and other velocity-dependent effects.
///
/// # Coordinate Convention
///
/// All velocities are returned in **ecliptic coordinates** (J2000.0 mean ecliptic)
/// with units of **AU per day**.
pub trait VelocityEphemeris: BodyEphemeris {
    /// Returns the velocity of `body` relative to the solar system barycenter.
    ///
    /// # Arguments
    ///
    /// - `body`: The body to query.
    /// - `jd`: The Julian Date (TDB scale preferred).
    ///
    /// # Returns
    ///
    /// Velocity `[vx, vy, vz]` in AU/day, in ecliptic J2000.0 coordinates.
    fn velocity_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3];

    /// Returns the velocity of `body` relative to the Sun (heliocentric).
    ///
    /// # Default Implementation
    ///
    /// Computes `velocity_barycentric(body) - velocity_barycentric(Sun)`.
    fn velocity_heliocentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        let [bvx, bvy, bvz] = self.velocity_barycentric(body, jd);
        let [svx, svy, svz] = self.velocity_barycentric(BodyId::Sun, jd);
        [bvx - svx, bvy - svy, bvz - svz]
    }
}

/// Blanket implementation for references to ephemeris providers.
impl<E: BodyEphemeris + ?Sized> BodyEphemeris for &E {
    fn position_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        (*self).position_barycentric(body, jd)
    }

    fn position_heliocentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        (*self).position_heliocentric(body, jd)
    }
}

impl<E: VelocityEphemeris + ?Sized> VelocityEphemeris for &E {
    fn velocity_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        (*self).velocity_barycentric(body, jd)
    }

    fn velocity_heliocentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        (*self).velocity_heliocentric(body, jd)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    /// A mock ephemeris that returns fixed positions for testing.
    ///
    /// This allows tests to verify that transforms use the ephemeris trait
    /// rather than hardcoded VSOP87 calls.
    #[derive(Debug, Clone)]
    pub struct MockEphemeris {
        pub sun_pos: [f64; 3],
        pub earth_pos: [f64; 3],
    }

    impl Default for MockEphemeris {
        fn default() -> Self {
            Self {
                sun_pos: [0.001, 0.002, 0.0001],
                earth_pos: [1.0, 0.0, 0.0],
            }
        }
    }

    impl BodyEphemeris for MockEphemeris {
        fn position_barycentric(&self, body: BodyId, _jd: JulianDate) -> [f64; 3] {
            match body {
                BodyId::Sun => self.sun_pos,
                BodyId::Earth => self.earth_pos,
                _ => [0.0, 0.0, 0.0],
            }
        }
    }

    #[test]
    fn test_heliocentric_from_barycentric() {
        let eph = MockEphemeris {
            sun_pos: [0.001, 0.002, 0.0001],
            earth_pos: [1.0, 0.0, 0.0],
        };

        let earth_helio = eph.position_heliocentric(BodyId::Earth, JulianDate::J2000);

        // heliocentric = barycentric - sun_barycentric
        assert!((earth_helio[0] - 0.999).abs() < 1e-10);
        assert!((earth_helio[1] - (-0.002)).abs() < 1e-10);
        assert!((earth_helio[2] - (-0.0001)).abs() < 1e-10);
    }

    #[test]
    fn test_body_id_display() {
        assert_eq!(format!("{}", BodyId::Sun), "Sun");
        assert_eq!(format!("{}", BodyId::Earth), "Earth");
    }
}
