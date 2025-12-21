//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! This module re-exports the algebraic types from `affn` and provides
//! frame-specific extension traits with convenient constructors and accessors
//! (e.g., `new_icrs(ra, dec)`, `ra()`, `dec()`).
//!
//! # Usage
//!
//! Import the extension traits to access frame-specific methods:
//!
//! ```rust
//! use siderust::coordinates::spherical::ext::*;
//! use siderust::coordinates::spherical::direction::ICRS;
//! use qtty::*;
//!
//! // Extension trait method for ICRS
//! let dir = ICRS::new_icrs(120.0 * DEG, 45.0 * DEG);
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! assert_eq!(dir.dec(), 45.0 * DEG);
//! ```

use crate::coordinates::{centers, frames};

pub use affn::spherical::{Direction, Position};

// =============================================================================
// Direction type aliases (frame-only, no center)
// =============================================================================

pub mod direction {
    pub use super::Direction;
    use super::frames;

    /// **Ecliptic** direction (longitude *L*, latitude *B*).
    pub type Ecliptic = Direction<frames::Ecliptic>;
    /// **Equatorial** direction (right‑ascension *α*, declination *δ*).
    pub type Equatorial = Direction<frames::Equatorial>;
    /// **Horizontal** direction (azimuth *Az*, altitude *Alt*).
    pub type Horizontal = Direction<frames::Horizontal>;
    /// **ICRS** direction.
    pub type ICRS = Direction<frames::ICRS>;
    /// **Geographic** (ECEF) direction: latitude, longitude.
    pub type Geographic = Direction<frames::ECEF>;
}

// =============================================================================
// Position type aliases (center + frame + distance)
// =============================================================================

pub mod position {
    pub use super::Position;
    use super::{centers, frames};
    use qtty::Kilometer;

    /// **Heliocentric Ecliptic** coordinates *(L, B, R)*.
    ///
    /// * `L` – ecliptic longitude, degrees in `[0, 360)`
    /// * `B` – ecliptic latitude,  degrees in `[-90, 90]`
    /// * `R` – heliocentric distance in unit `U` (e.g. `AstronomicalUnit`)
    pub type Ecliptic<U> = Position<centers::Heliocentric, frames::Ecliptic, U>;

    /// **Geocentric Equatorial** coordinates *(δ, α, d)*.
    ///
    /// * `δ` – declination, degrees in `[-90, 90]`
    /// * `α` – right‑ascension, degrees in `[0, 360)`
    /// * `d` – geocentric distance in unit `U` (e.g. `Kilometer`)
    pub type Equatorial<U> = Position<centers::Geocentric, frames::Equatorial, U>;

    /// **Topocentric Horizontal** coordinates *(Alt, Az, d)*.
    ///
    /// * `Alt` – altitude above the horizon, degrees in `[-90, 90]`
    /// * `Az`  – azimuth from the north, degrees in `[0, 360)`
    /// * `d`   – straight‑line distance from the observer in unit `U`
    pub type Horizontal<U> = Position<centers::Topocentric, frames::Horizontal, U>;

    /// **Barycentric ICRS** coordinates.
    pub type ICRS<U> = Position<centers::Barycentric, frames::ICRS, U>;
    /// **Heliocentric ICRS** coordinates.
    pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
    /// **Geocentric ICRS** coordinates.
    pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;

    /// **Geographic (ECEF)** position: latitude, longitude, altitude (km).
    pub type Geographic = Position<centers::Geocentric, frames::ECEF, Kilometer>;
}
