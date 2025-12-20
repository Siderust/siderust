//! # Position‐type Specialisations
//!
//! This module defines the core spherical [`Position`] type (center + frame + distance)
//! and provides domain‑specific type aliases
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

use crate::coordinates::algebra::{centers, frames};
use qtty::*;

use std::marker::PhantomData;

/// A spherical **position** (center + frame + distance).
///
/// This is the fundamental spherical coordinate type used across the crate.
/// Spherical directions are represented separately by [`super::direction::Direction<F>`]
/// and intentionally have **no** reference center.
#[derive(Debug, Clone, Copy)]
pub struct Position<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: Unit> {
    pub polar: Degrees,        // θ (polar/latitude/declination)
    pub azimuth: Degrees,      // φ (azimuth/longitude/right ascension)
    pub distance: Quantity<U>, // radial distance

    center_params: C::Params,
    _frame: PhantomData<F>,
}

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Constructs a new spherical position with explicit center parameters.
    pub const fn new_raw_with_params(
        center_params: C::Params,
        polar: Degrees,
        azimuth: Degrees,
        distance: Quantity<U>,
    ) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            center_params,
            _frame: PhantomData,
        }
    }

    /// Returns a reference to the center parameters.
    pub fn center_params(&self) -> &C::Params {
        &self.center_params
    }

    /// Calculates the angular separation between this position and another.
    pub fn angular_separation(&self, other: Self) -> Degrees {
        let az1 = self.azimuth.to::<Radian>();
        let po1 = self.polar.to::<Radian>();
        let az2 = other.azimuth.to::<Radian>();
        let po2 = other.polar.to::<Radian>();

        let x = (po1.cos() * po2.sin()) - (po1.sin() * po2.cos() * (az2 - az1).cos());
        let y = po2.cos() * (az2 - az1).sin();
        let z = (po1.sin() * po2.sin()) + (po1.cos() * po2.cos() * (az2 - az1).cos());

        let angle_rad = (x * x + y * y).sqrt().atan2(z);
        Radians::new(angle_rad).to::<Degree>()
    }

    /// Extracts the corresponding spherical **direction** (frame-only).
    #[must_use]
    pub fn direction(&self) -> super::direction::Direction<F> {
        super::direction::Direction::new(self.polar, self.azimuth)
    }
}

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Convenience constructor for centers with `Params = ()`.
    pub const fn new_raw(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self {
        Self::new_raw_with_params((), polar, azimuth, distance)
    }
}

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
        // new_raw(polar, azimuth, distance)
        let coord = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(Degrees::new(90.0), Degrees::new(45.0), 1.0 * AU);
        assert_eq!(coord.polar.value(), 90.0);
        assert_eq!(coord.azimuth.value(), 45.0);
        assert_eq!(coord.distance.value(), 1.0);
    }

    #[test]
    fn test_spherical_coord_to_string() {
        // new_raw(polar, azimuth, distance)
        let coord = Position::<centers::Geocentric, frames::ICRS, AstronomicalUnit>::new_raw(Degrees::new(60.0), Degrees::new(30.0), 1000.0 * AU);
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 60"));
        assert!(coord_string.contains("φ: 30"));
        assert!(coord_string.contains("r: 1000"));
    }

    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = Position::<centers::Heliocentric, frames::ICRS, AstronomicalUnit>::new_raw(0.0 * DEG, 0.0 * DEG, 0.0 * AU);
        assert_eq!(coord.polar.value(), 0.0);
        assert_eq!(coord.azimuth.value(), 0.0);
        assert_eq!(coord.distance.value(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        // new_raw(polar, azimuth, distance)
        let coord = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(45.123456 * DEG, 90.654321 * DEG, 1234.56789 * AU);
        assert!((coord.polar.value() - 45.123456).abs() < 1e-6);
        assert!((coord.azimuth.value() - 90.654321).abs() < 1e-6);
        assert!((coord.distance - 1234.56789 * AU).abs() < 1e-6 * AU);
    }

    #[test]
    fn direction_returns_unit_vector() {
        let pos = Position::<centers::Heliocentric, frames::Ecliptic, AstronomicalUnit>::new_raw(10.0 * DEG, 20.0 * DEG, 2.5 * AU);
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

        let c = Position::<centers::Geocentric, frames::Equatorial, Kilometer>::CENTER;
        assert_eq!(c.polar.value(), 0.0);
        assert_eq!(c.azimuth.value(), 0.0);
        assert_eq!(c.distance.value(), 0.0);
    }

    #[test]
    fn from_degrees_matches_new_raw() {
        let a = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(45.0 * DEG, 30.0 * DEG, 3.0 * AU);
        let b = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(45.0 * DEG, 30.0 * DEG, 3.0 * AU);
        assert_eq!(a.polar, b.polar);
        assert_eq!(a.azimuth, b.azimuth);
        assert_eq!(a.distance, b.distance);
    }

    #[test]
    fn distance_identity_zero_and_orthogonal() {
        // identity
        let a = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(0.0 * DEG, 0.0 * DEG, 1.0 * AU);
        let d0 = a.distance_to(&a);
        assert!(d0.abs().value() < EPS);

        // orthogonal points on unit sphere → chord length sqrt(2) * r
        let b = Position::<centers::Barycentric, frames::ICRS, AstronomicalUnit>::new_raw(0.0 * DEG, 90.0 * DEG, 1.0 * AU);
        let d = a.distance_to(&b);
        assert!((d.value() - SQRT_2).abs() < EPS);
    }
}
