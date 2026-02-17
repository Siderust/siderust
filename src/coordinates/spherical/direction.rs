//! Frame-specific direction type aliases.
//!
//! These provide convenient shorthand for common direction types.

use super::Direction;
use crate::coordinates::frames;

/// **EclipticMeanJ2000** direction (longitude *λ*, latitude *β*).
pub type EclipticMeanJ2000 = Direction<frames::EclipticMeanJ2000>;
/// **Equatorial mean J2000** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialMeanJ2000 = Direction<frames::EquatorialMeanJ2000>;
/// **Equatorial mean of date** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialMeanOfDate = Direction<frames::EquatorialMeanOfDate>;
/// **Equatorial true of date** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialTrueOfDate = Direction<frames::EquatorialTrueOfDate>;
/// **Horizontal** direction (altitude *Alt*, azimuth *Az*).
///
/// Azimuth follows the **North-clockwise** convention (0° = North, increasing
/// through East). For data that uses a different convention, see
/// [`crate::coordinates::horizontal`] for conversion helpers.
pub type Horizontal = Direction<frames::Horizontal>;
/// **ICRS** direction.
pub type ICRS = Direction<frames::ICRS>;
/// **Geographic** (ECEF) direction: longitude, latitude.
pub type Geographic = Direction<frames::ECEF>;
/// **Galactic** direction (l, b).
pub type Galactic = Direction<frames::Galactic>;
