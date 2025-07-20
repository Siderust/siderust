//! # Spherical Coordinates
//!
//! This module defines the generic [`SphericalCoord<C, F>`] type for representing positions or directions
//! in spherical coordinates, parameterized by astronomical reference centers and frames for strong type safety.
//!
//! ## Overview
//!
//! - **Generic over Center and Frame:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`, `Barycentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Angles are stored as [`Degrees`]; distance is optional and typically in AstronomicalUnits or parsecs (see context).
//! - **Conversions:** Seamless conversion to and from [`Vector`] via `From`/`Into`.
//! - **Operations:** Compute Euclidean distance and angular separation between coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::{Degrees, AstronomicalUnit};
//!
//! // Create a heliocentric ecliptic spherical position
//! let sph = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(Degrees::new(45.0), Degrees::new(7.0), 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance).
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use crate::units::*;
use crate::coordinates::{frames, centers};

use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`).
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`).
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
> {
    pub polar: Degrees,        // θ (polar/latitude/declination)
    pub azimuth: Degrees,      // φ (azimuth/longitude/right ascension)
    pub distance: Quantity<U>, // Distance (AstronomicalUnits, parsec, etc.)

    _center: PhantomData<C>,
    _frame: PhantomData<F>,
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Constructs a new spherical coordinate.
    ///
    /// * `polar`: angle from the reference plane, in degrees  
    /// * `azimuth`: angle from the reference direction, in degrees  
    /// * `distance`: radial distance in the same unit `U`
    pub const fn new_raw(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            _center: PhantomData,
            _frame: PhantomData,
        }
    }

    /// Calculates the angular separation between this coordinate and another.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The angular separation in degrees.
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
}

#[cfg(test)]
mod tests {
    use crate::coordinates::spherical::SphericalCoord;
    use crate::coordinates::centers::Geocentric;
    use crate::coordinates::frames::ICRS;
    use crate::units::*;

    const EPS: Degrees = Degrees::new(1e-6);       // tolerance for the exact geometry cases
    const EPS_STAR: Degrees = Degrees::new(1e-2);  // tolerance for the real‑world catalogue case

    /// Helper to build a unit‑length direction in the ICRS frame.
    fn dir(dec: f64, ra: f64) -> SphericalCoord<Geocentric, ICRS, Unitless> {
        SphericalCoord::new_raw(
            Degrees::new(dec),   // polar / declination
            Degrees::new(ra),    // azimuth / right‑ascension
            Quantity::<Unitless>::new(1.0),
        )
    }

    #[test]
    fn identity_separation_is_zero() {
        let a = dir(12.3, 45.6);
        let sep = a.angular_separation(a);
        assert!(sep.to::<Degree>().abs() < EPS, "expected 0°, got {}°", sep);
    }

    #[test]
    fn orthogonal_points_give_ninety_degrees() {
        let a = dir(0.0, 0.0);
        let b = dir(0.0, 90.0);
        let sep = a.angular_separation(b);
        assert!((sep.to::<Degree>() - 90.0*DEG).abs() < EPS, "expected 90°, got {}°", sep);
    }

    #[test]
    fn antipodal_points_give_180_degrees() {
        let a = dir(0.0, 0.0);
        let b = dir(0.0, 180.0);
        let sep = a.angular_separation(b);
        assert!((sep.to::<Degree>() - 180.0*DEG).abs() < EPS, "expected 180°, got {}°", sep);
    }

    #[test]
    fn polaris_betelgeuse_real_world() {
        // Star coordinates (J2000) from SIMBAD
        let polaris = dir(89.26410897, 37.95456067);   // Dec, RA
        let betel   = dir(7.407064,    88.792939);     // Dec, RA
        let sep = polaris.angular_separation(betel);
        assert!((sep.to::<Degree>() - 82.1286*DEG).abs() < EPS_STAR,
            "expected 224882.13°, got {}°", sep);
    }

}