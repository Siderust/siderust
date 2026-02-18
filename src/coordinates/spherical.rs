// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! The core spherical coordinate functionality lives in the `affn` crate
//! (domain-agnostic geometry kernel). This module provides:
//!
//! - Astronomical reference frame types (ICRS, EclipticMeanJ2000, Horizontal, etc.)
//! - Inherent named constructors and getters (`ra`/`dec`, `lon`/`lat`, `alt`/`az`)
//! - Convenient type aliases for common coordinate combinations

// Re-export affn types directly — no wrappers
pub use affn::spherical::Direction;
pub use affn::spherical::Position;

/// Direction aliases are grouped under the `direction` submodule to avoid
/// name collisions with position aliases (which live in `position`).
pub mod direction {
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
	/// **ECEF** direction: unit vector in the Earth-fixed frame.
	///
	/// For geodetic (lon/lat/h) positions, use [`Geodetic::<ECEF>`](crate::coordinates::centers::Geodetic) instead;
	/// this type is for unit vectors, not geodetic positions.
	pub type EcefDir = Direction<frames::ECEF>;
	/// **Galactic** direction (l, b).
	pub type Galactic = Direction<frames::Galactic>;
}

/// Position aliases are grouped under the `position` submodule.
pub mod position {
	//! Frame and center-specific position type aliases.
	//!
	//! These provide convenient shorthand for common position types.

	use super::Position;
	use crate::coordinates::{centers, frames};

	/// **Heliocentric EclipticMeanJ2000** coordinates *(λ, β, R)*.
	///
	/// * `λ` – ecliptic longitude, degrees in `[0, 360)`
	/// * `β` – ecliptic latitude,  degrees in `[-90, 90]`
	/// * `R` – heliocentric distance in unit `U` (e.g. `AstronomicalUnit`)
	pub type EclipticMeanJ2000<U> = Position<centers::Heliocentric, frames::EclipticMeanJ2000, U>;

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
	///
	/// Azimuth follows the **North-clockwise** convention (0° = North, increasing
	/// through East). For data that uses a different convention, see
	/// [`crate::coordinates::horizontal`] for conversion helpers.
	pub type Horizontal<U> = Position<centers::Topocentric, frames::Horizontal, U>;

	/// **Barycentric ICRS** coordinates.
	pub type ICRS<U> = Position<centers::Barycentric, frames::ICRS, U>;
	/// **Heliocentric ICRS** coordinates.
	pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
	/// **Geocentric ICRS** coordinates.
	///
	/// # Approximation
	///
	/// This alias uses [`frames::ICRS`] as a first-order approximation for
	/// the Geocentric Celestial Reference System ([`frames::GCRS`]). The
	/// difference is < 1 mas for typical astronomical directions (neglected:
	/// geocentre offset, relativistic terms). For strictly IAU-correct GCRS,
	/// use `Position<Geocentric, frames::GCRS, U>` directly.
	pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;

	// NOTE: The `Geographic` spherical-position alias has been removed.
	// The old definition `Geographic = Position<Geocentric, ECEF, Kilometer>`
	// was a correctness footgun: in a spherical position `distance` is radial
	// distance, not ellipsoidal height, so `.to_cartesian()` produced
	// geometrically wrong results for geodetic coordinates.
	//
	// Use `affn::ellipsoidal::Position` for geodetic (lon/lat/h) constants
	// and `to_cartesian()` for the ellipsoid-correct geodetic -> ECEF conversion.
}
