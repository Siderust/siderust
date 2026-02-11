//! Frame-specific direction type aliases.
//!
//! These provide convenient shorthand for common direction types.

use super::Direction;
use crate::coordinates::frames;

/// **Ecliptic** direction (longitude *λ*, latitude *β*).
pub type Ecliptic = Direction<frames::Ecliptic>;
/// **Equatorial mean J2000** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialMeanJ2000 = Direction<frames::EquatorialMeanJ2000>;
/// **Equatorial mean of date** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialMeanOfDate = Direction<frames::EquatorialMeanOfDate>;
/// **Equatorial true of date** direction (right‑ascension *α*, declination *δ*).
pub type EquatorialTrueOfDate = Direction<frames::EquatorialTrueOfDate>;
/// **Horizontal** direction (altitude *Alt*, azimuth *Az*).
pub type Horizontal = Direction<frames::Horizontal>;
/// **ICRS** direction.
pub type ICRS = Direction<frames::ICRS>;
/// **Geographic** (ECEF) direction: longitude, latitude.
pub type Geographic = Direction<frames::ECEF>;
/// **Galactic** direction (l, b).
pub type Galactic = Direction<frames::Galactic>;
