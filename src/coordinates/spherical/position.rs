//! # Position‐type Specialisations
//!
//! This module re‑exports [`SphericalCoord`] under domain‑specific type aliases
//! (e.g. *heliocentric ecliptic position* or *geocentric equatorial position*)
//! and provides convenience constructors and helpers that apply equally to
//! every centre / frame combination.
//!
//! Each alias fixes the reference **centre** and **frame** at the type level
//! while leaving the *distance unit* `U` generic, so you can work in
//! `AstronomicalUnit`, `Kilometer`, or any other [`LengthUnit`] supported by
//! `siderust::units`.
//!
//! ## Quick example
//! ```rust
//! use siderust::coordinates::spherical::position::{Ecliptic, Equatorial};
//! use qtty::*;
//!
//! // Heliocentric ecliptic coordinates of Earth at J2000
//! let ecl = Ecliptic::<AstronomicalUnit>::new(0.0, 0.0, 1.0*AU);
//!
//! // Convert to a direction (unit vector)
//! let dir = ecl.direction();
//! println!("Sun–Earth direction = {dir}");
//! ```
//!
//! The aliases are:
//! | Alias | Centre | Frame | Components (θ, φ, r) |
//! |-------|--------|-------|-----------------------|
//! | `Ecliptic<U>`   | `Heliocentric` | `Ecliptic`   | longitude *L*, latitude *B*, radius *R* |
//! | `Equatorial<U>` | `Geocentric`   | `Equatorial` | declination *δ*, right‑ascension *α*, distance *d* |
//! | `Horizontal<U>` | `Topocentric`  | `Horizontal` | altitude *Alt*, azimuth *Az*, distance *d* |
//! | `ICRS<U>`       | `Barycentric`  | `ICRS`       | declination *δ*, right‑ascension *α*, distance *d* |
//! | `HCRS<U>`       | `Heliocentric` | `ICRS`       | declination *δ*, right‑ascension *α*, distance *d* |
//! | `GCRS<U>`       | `Geocentric`   | `ICRS`       | declination *δ*, right‑ascension *α*, distance *d* |
//! | `Geographic`    | `Geocentric`   | `ECEF`       | latitude *φ*, longitude *λ*, altitude *h* (km) |
//!
//! All aliases are **zero‑cost** at compile time: they are simple `type` synonyms.

use super::SphericalCoord;
use crate::coordinates::{centers, frames};
use qtty::*;

// TODO: Bound U to LengthUnit and VelocityUnit
// see issue #112792 <https://github.com/rust-lang/rust/issues/112792> for more information
//pub type Position<C, F, U: LengthUnit>  = SphericalCoord<C, F, U>;

/// Generic *position* alias: a spherical coordinate that **does** carry radial
/// distance information.
///
/// You will seldom use this directly; reach for a more specific alias such as
/// [`Ecliptic`], [`Equatorial`], etc., to get compile‑time guarantees about the
/// reference centre / frame.
pub type Position<C, F, U> = SphericalCoord<C, F, U>;

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

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{
    /// The *origin* of this coordinate system (all angles 0, radius 0). AKA Null Vector.
    pub const CENTER: Self = Self::new_raw(
        Degrees::new(0.0),
        Degrees::new(0.0),
        Quantity::<U>::new(0.0),
    );
}

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{

    /// Euclidean distance to another position **in the same centre & frame**.
    ///
    /// The result is expressed in the *same unit `U`* as the inputs.
    #[must_use]
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        self.to_cartesian().distance_to(&other.to_cartesian())
    }
}

impl<C, F, U> std::fmt::Display for Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}, r: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.polar,
            self.azimuth,
            self.distance
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::SQRT_2;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_spherical_coord_creation() {
        let coord = ICRS::<AstronomicalUnit>::new(Degrees::new(45.0), Degrees::new(90.0), 1.0);
        assert_eq!(coord.ra().value(), 45.0);
        assert_eq!(coord.dec().value(), 90.0);
        assert_eq!(coord.distance.value(), 1.0);
    }

    #[test]
    fn test_spherical_coord_to_string() {
        let coord = GCRS::<AstronomicalUnit>::new(Degrees::new(30.0), Degrees::new(60.0), 1000.0);
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 60"));
        assert!(coord_string.contains("φ: 30"));
        assert!(coord_string.contains("r: 1000"));
    }

    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = HCRS::<AstronomicalUnit>::new(0.0 * DEG, 0.0 * DEG, 0.0);
        assert_eq!(coord.polar.value(), 0.0);
        assert_eq!(coord.azimuth.value(), 0.0);
        assert_eq!(coord.distance.value(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = ICRS::<AstronomicalUnit>::new(90.654321 * DEG, 45.123456 * DEG, 1234.56789);
        assert!((coord.dec().value() - 45.123456).abs() < 1e-6);
        assert!((coord.ra().value() - 90.654321).abs() < 1e-6);
        assert!((coord.distance - 1234.56789 * AU).abs() < 1e-6 * AU);
    }

    #[test]
    fn direction_returns_unit_vector() {
        let pos = Ecliptic::<AstronomicalUnit>::new_raw(10.0 * DEG, 20.0 * DEG, 2.5 * AU);
        let dir = pos.direction();

        // radial component must be exactly 1 (unitless)
        assert_eq!(dir.distance.value(), 1.0);

        // angular components are preserved
        assert!((dir.polar - 10.0 * DEG).abs() < EPS * DEG);
        assert!((dir.azimuth - 20.0 * DEG).abs() < EPS * DEG);
    }

    #[test]
    fn center_constant_is_origin() {
        use qtty::Kilometer;

        let c = Equatorial::<Kilometer>::CENTER;
        assert_eq!(c.polar.value(), 0.0);
        assert_eq!(c.azimuth.value(), 0.0);
        assert_eq!(c.distance.value(), 0.0);
    }

    #[test]
    fn from_degrees_matches_new_raw() {
        let a = ICRS::<AstronomicalUnit>::new(45.0 * DEG, 30.0 * DEG, 3.0 * AU);
        let b = ICRS::<AstronomicalUnit>::new(45.0 * DEG, 30.0 * DEG, 3.0 * AU);
        assert_eq!(a.polar, b.polar);
        assert_eq!(a.azimuth, b.azimuth);
        assert_eq!(a.distance, b.distance);
    }

    #[test]
    fn distance_identity_zero_and_orthogonal() {
        // identity
        let a = ICRS::<AstronomicalUnit>::new(0.0 * DEG, 0.0 * DEG, 1.0 * AU);
        let d0 = a.distance_to(&a);
        assert!(d0.abs().value() < EPS);

        // orthogonal points on unit sphere → chord length sqrt(2) * r
        let b = ICRS::<AstronomicalUnit>::new(0.0 * DEG, 90.0 * DEG, 1.0 * AU);
        let d = a.distance_to(&b);
        assert!((d.value() - SQRT_2).abs() < EPS);
    }
}
