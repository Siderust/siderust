//! # Direction‐type Specialisations (unit vectors)
//!
//! A **direction** represents a *unit‐length* pointing vector: spherical
//! coordinates in which the radial component is implicitly fixed to 1.
//! Instead of carrying a distance, the type alias [`Direction<C, F>`] hard‐codes
//! the distance unit to [`Unitless`].  This allows you to express pointing
//! vectors (e.g. *line of sight*, *observer zenith*, etc.) while still keeping
//! the **centre** and **frame** strongly typed.
//!
//! Most users will prefer the ready‑made aliases below, each corresponding to a
//! widely used astronomical system:
//!
//! | Alias        | Centre           | Frame        | Angles (θ, φ)                        |
//! |--------------|------------------|--------------|--------------------------------------|
//! | `Ecliptic`   | `Heliocentric`   | `Ecliptic`   | latitude *B*, longitude *L*          |
//! | `Equatorial` | `Geocentric`     | `Equatorial` | declination *δ*, right‑ascension *α* |
//! | `Horizontal` | `Topocentric`    | `Horizontal` | altitude *Alt*, azimuth *Az*         |
//! | `ICRS`       | `Barycentric`    | `ICRS`       | declination *δ*, right‑ascension *α* |
//! | `HCRS`       | `Heliocentric`   | `ICRS`       | declination *δ*, right‑ascension *α* |
//! | `GCRS`       | `Geocentric`     | `ICRS`       | declination *δ*, right‑ascension *α* |
//! | `Geographic` | `Geocentric`     | `ECEF`       | latitude *φ*, longitude *λ*          |
//!
//! ## Quick example
//! ```rust
//! use siderust::coordinates::spherical::direction::ICRS;
//! use qtty::*;
//!
//! // Barycentric pointing direction to Vega (J2000)
//! let vega: ICRS = ICRS::new(279.23473479 * DEG, 38.78368896 * DEG);
//!
//! // Convert that direction into a position one parsec away:
//! use qtty::Parsec;
//! let one_pc = 1.0 * PC;
//! let pos = vega.position(one_pc);
//! println!("Vega direction  = {vega}\nVega@1pc position = {pos}");
//! ```
//!
//! The conversion above is *zero‑cost*; both `Direction` and `Position` are
//! simple aliases around [`SphericalCoord`]. Only the *type* (and therefore the
//! allowed operations) change at compile time — there is no run‑time penalty.
//!
//! ---
//!
//! ## API summary
//! * `position(magnitude)` – promotes a direction to a position with the given
//!   radial `magnitude`.
//! * `Display` – prints the centre, frame and angles.

use super::SphericalCoord;
use crate::coordinates::{centers, centers::ReferenceCenter, frames, frames::ReferenceFrame};
use qtty::{Dimension, LengthUnit, Quantity, Unit};

/// Marker dimension for direction (dimensionless unit vector).
pub enum DirectionDim {}
impl Dimension for DirectionDim {}

/// Marker unit for direction types (unit vectors with implicit radius = 1).
///
/// This is a local type distinct from qtty's `Unitless`, which allows
/// Direction and Position to have non-overlapping impl blocks.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct DirectionUnit;

impl Unit for DirectionUnit {
    const RATIO: f64 = 1.0;
    type Dim = DirectionDim;
    const SYMBOL: &'static str = "";
}

/// Generic alias for a *unit vector* (radius = 1).
///
/// The distance unit is fixed to [`DirectionUnit`].
pub type Direction<C, F> = SphericalCoord<C, F, DirectionUnit>;

/// **Heliocentric ecliptic** direction (longitude *L*, latitude *B*).
pub type Ecliptic = Direction<centers::Heliocentric, frames::Ecliptic>;
/// **Geocentric equatorial** direction (right‑ascension *α*, declination *δ*).
pub type Equatorial = Direction<centers::Geocentric, frames::Equatorial>;
/// **Topocentric horizontal** direction (azimuth *Az*, altitude *Alt*).
pub type Horizontal = Direction<centers::Topocentric, frames::Horizontal>;
/// **Barycentric ICRS** direction.
pub type ICRS = Direction<centers::Barycentric, frames::ICRS>;
/// **Heliocentric ICRS** direction.
pub type HCRS = Direction<centers::Heliocentric, frames::ICRS>;
/// **Geocentric ICRS** direction.
pub type GCRS = Direction<centers::Geocentric, frames::ICRS>;
/// **Geographic** (ECEF) direction: latitude, longitude.
pub type Geographic = Direction<centers::Geocentric, frames::ECEF>;

impl<C, F> Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    /// Promotes this direction to a full [`Position`](super::Position) with the
    /// supplied radial `magnitude`.
    ///
    /// # Example
    /// ```rust
    /// use siderust::coordinates::spherical::direction::Ecliptic;
    /// use qtty::*;
    ///
    /// let dir = Ecliptic::new(0.0*DEG, 0.0*DEG);
    /// let pos = dir.position(1.0*AU); // 1 au
    /// assert_eq!(pos.distance.value(), 1.0);
    /// ```
    #[must_use]
    pub fn position<U: LengthUnit>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U> {
        super::Position::new_raw(self.polar, self.azimuth, magnitude)
    }
}

impl<C, F> std::fmt::Display for Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.polar,
            self.azimuth
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn creates_valid_spherical_position() {
        let polar = Degrees::new(45.0);
        let azimuth = Degrees::new(90.0);

        let coord = ICRS::new(polar, azimuth);

        assert_eq!(coord.ra().value(), 45.0);
        assert_eq!(coord.dec().value(), 90.0);
        assert_eq!(coord.distance.value(), 1.0);
    }

    #[test]
    fn displays_coordinate_as_string_correctly() {
        let coord = ICRS::new(Degrees::new(30.0), Degrees::new(60.0));

        let output = coord.to_string();

        assert!(output.contains("θ: 60"), "Missing polar angle");
        assert!(output.contains("φ: 30"), "Missing azimuth");
    }

    #[test]
    fn maintains_high_precision_on_values() {
        let polar = Degrees::new(90.654_321);
        let azimuth = Degrees::new(45.123_456);

        let coord = ICRS::new(polar, azimuth);

        assert!((coord.ra().value() - 90.654_321).abs() < 1e-6);
        assert!((coord.dec().value() - 45.123_456).abs() < 1e-6);
    }

    const EPS: f64 = 1e-6;

    #[test]
    fn position_method_promotes_with_given_radius() {
        let dir: ICRS = ICRS::new(Degrees::new(120.0), Degrees::new(-30.0));
        let pos = dir.position(Quantity::<AstronomicalUnit>::new(2.0));

        // angles are preserved
        assert!((pos.ra().value() - 120.0).abs() < EPS);
        assert!((pos.dec().value() - (-30.0)).abs() < EPS);

        // distance matches the supplied magnitude
        assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
    }

    #[test]
    fn direction_display_mentions_center_and_frame() {
        let eq: Equatorial = Equatorial::new(Degrees::new(45.0), Degrees::new(10.0));
        let s = eq.to_string();
        assert!(s.contains("Geocentric"), "missing center");
        assert!(s.contains("Equatorial"), "missing frame");
    }
}
