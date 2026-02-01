//! Frame and center-specific position type aliases.
//!
//! These provide convenient shorthand for common position types.

use super::Position;
use crate::coordinates::{centers, frames};
use qtty::Kilometer;

/// **Heliocentric Ecliptic** coordinates *(λ, β, R)*.
///
/// * `λ` – ecliptic longitude, degrees in `[0, 360)`
/// * `β` – ecliptic latitude,  degrees in `[-90, 90]`
/// * `R` – heliocentric distance in unit `U` (e.g. `AstronomicalUnit`)
pub type Ecliptic<U> = Position<centers::Heliocentric, frames::Ecliptic, U>;

/// **Geocentric Equatorial mean J2000** coordinates *(α, δ, d)*.
///
/// * `α` – right‑ascension, degrees in `[0, 360)`
/// * `δ` – declination, degrees in `[-90, 90]`
/// * `d` – geocentric distance in unit `U` (e.g. `Kilometer`)
pub type EquatorialMeanJ2000<U> = Position<centers::Geocentric, frames::EquatorialMeanJ2000, U>;

/// **Geocentric Equatorial mean of date** coordinates *(α, δ, d)*.
pub type EquatorialMeanOfDate<U> = Position<centers::Geocentric, frames::EquatorialMeanOfDate, U>;

/// **Geocentric Equatorial true of date** coordinates *(α, δ, d)*.
pub type EquatorialTrueOfDate<U> = Position<centers::Geocentric, frames::EquatorialTrueOfDate, U>;

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

/// **Geographic (ECEF)** position: longitude, latitude, altitude (km).
pub type Geographic = Position<centers::Geocentric, frames::ECEF, Kilometer>;
