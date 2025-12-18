//! # Direction‐type Specialisations (unit vectors)
//!
//! A **direction** represents a *unit‐length* pointing vector: spherical
//! coordinates in which the radial component is implicitly fixed to 1.
//!
//! ## Mathematical Foundations
//!
//! Directions are **free vectors**: they live in the vector space, not in affine
//! space. Unlike positions, directions are translation-invariant:
//!
//! - Rotating a direction is valid (frame transformation)
//! - "Translating" a direction to a different origin is **meaningless**
//!
//! Therefore, `Direction<F>` has no center parameter, only a frame `F`.
//!
//! ## Observer-Dependent Directions
//!
//! To compute the direction from an observer to a target (line of sight), use
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) with two
//! positions. This correctly accounts for the observer's location.
//!
//! ## Type Aliases
//!
//! | Alias        | Frame        | Angles (θ, φ)                        |
//! |--------------|--------------|--------------------------------------|
//! | `Ecliptic`   | `Ecliptic`   | latitude *B*, longitude *L*          |
//! | `Equatorial` | `Equatorial` | declination *δ*, right‑ascension *α* |
//! | `Horizontal` | `Horizontal` | altitude *Alt*, azimuth *Az*         |
//! | `ICRS`       | `ICRS`       | declination *δ*, right‑ascension *α* |
//! | `Geographic` | `ECEF`       | latitude *φ*, longitude *λ*          |
//!
//! ## Quick example
//! ```rust
//! use siderust::coordinates::spherical::direction::ICRS;
//! use qtty::*;
//!
//! // Direction to Vega (J2000)
//! let vega: ICRS = ICRS::new(279.23473479 * DEG, 38.78368896 * DEG);
//!
//! // Convert that direction into a position one parsec away:
//! use siderust::coordinates::centers::Barycentric;
//! let one_pc = 1.0 * PC;
//! let pos = vega.position::<Barycentric, _>(one_pc);
//! println!("Vega direction = {vega}\nVega@1pc position = {pos}");
//! ```

use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::{self, ReferenceFrame};
use qtty::{Degrees, Dimension, LengthUnit, Quantity, Unit};

use std::marker::PhantomData;

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

/// A spherical direction (unit vector) in a specific reference frame.
///
/// Directions are frame-dependent but center-independent (free vectors).
/// They cannot undergo center transformations, only frame transformations.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
#[derive(Debug, Clone, Copy)]
pub struct Direction<F: ReferenceFrame> {
    /// Polar angle (θ) - latitude, declination, or altitude, in degrees.
    pub polar: Degrees,
    /// Azimuthal angle (φ) - longitude, right ascension, or azimuth, in degrees.
    pub azimuth: Degrees,
    /// Distance (always 1 for unit vectors).
    pub distance: Quantity<DirectionUnit>,
    _frame: PhantomData<F>,
}

impl<F: ReferenceFrame> Direction<F> {
    /// Creates a new direction from polar and azimuth angles.
    ///
    /// * `polar`: angle from the reference plane, in degrees
    /// * `azimuth`: angle from the reference direction, in degrees
    pub fn new(polar: Degrees, azimuth: Degrees) -> Self {
        Self {
            polar,
            azimuth,
            distance: Quantity::<DirectionUnit>::new(1.0),
            _frame: PhantomData,
        }
    }

    /// Creates a direction from raw angle values.
    pub const fn new_raw(polar: Degrees, azimuth: Degrees) -> Self {
        Self {
            polar,
            azimuth,
            distance: Quantity::<DirectionUnit>::new(1.0),
            _frame: PhantomData,
        }
    }

    /// Promotes this direction to a full Position with the supplied radial magnitude.
    ///
    /// # Type Parameters
    /// - `C`: The center for the resulting position
    /// - `U`: The length unit for the magnitude
    ///
    /// # Example
    /// ```rust
    /// use siderust::coordinates::spherical::direction::Ecliptic;
    /// use siderust::coordinates::centers::Heliocentric;
    /// use qtty::*;
    ///
    /// let dir = Ecliptic::new(0.0*DEG, 0.0*DEG);
    /// let pos = dir.position::<Heliocentric, AstronomicalUnit>(1.0*AU);
    /// assert_eq!(pos.distance.value(), 1.0);
    /// ```
    #[must_use]
    pub fn position<C, U>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U>
    where
        C: ReferenceCenter<Params = ()>,
        U: LengthUnit,
    {
        super::Position::new_raw(self.polar, self.azimuth, magnitude)
    }

    /// Promotes this direction to a Position with explicit center parameters.
    #[must_use]
    pub fn position_with_params<C, U>(
        &self,
        center_params: C::Params,
        magnitude: Quantity<U>,
    ) -> super::Position<C, F, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
    {
        super::Position::new_raw_with_params(center_params, self.polar, self.azimuth, magnitude)
    }

    /// Returns the right ascension (alias for azimuth in RA/Dec systems).
    /// In spherical coordinates, RA is the azimuthal angle.
    pub fn ra(&self) -> Degrees {
        self.azimuth
    }

    /// Returns the declination (alias for polar in RA/Dec systems).
    /// In spherical coordinates, Dec is the polar angle.
    pub fn dec(&self) -> Degrees {
        self.polar
    }

    /// Calculates the angular separation between this direction and another.
    pub fn angular_separation(&self, other: &Self) -> Degrees {
        use qtty::{Radian, Radians, Degree};

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

    /// Converts to cartesian Direction.
    pub fn to_cartesian(&self) -> crate::coordinates::cartesian::Direction<F>
    where
        F: frames::MutableFrame,
    {
        use qtty::Radian;

        let polar_rad = self.polar.to::<Radian>();
        let azimuth_rad = self.azimuth.to::<Radian>();

        let x = azimuth_rad.cos() * polar_rad.cos();
        let y = azimuth_rad.sin() * polar_rad.cos();
        let z = polar_rad.sin();

        crate::coordinates::cartesian::Direction::<F>::new(x, y, z)
    }
}

impl<F: ReferenceFrame> std::fmt::Display for Direction<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Frame: {}, θ: {:.6}, φ: {:.6}",
            F::frame_name(),
            self.polar,
            self.azimuth
        )
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn creates_valid_spherical_direction() {
        let polar = Degrees::new(45.0);  // This is Dec
        let azimuth = Degrees::new(90.0); // This is RA

        let coord = ICRS::new(polar, azimuth);

        // In spherical coords: polar=Dec, azimuth=RA
        assert_eq!(coord.ra().value(), 90.0);  // RA = azimuth
        assert_eq!(coord.dec().value(), 45.0); // Dec = polar
        assert_eq!(coord.distance.value(), 1.0);
    }

    #[test]
    fn displays_coordinate_as_string_correctly() {
        let coord = ICRS::new(Degrees::new(30.0), Degrees::new(60.0));
        let output = coord.to_string();
        assert!(output.contains("θ: 30"), "Missing polar angle");
        assert!(output.contains("φ: 60"), "Missing azimuth");
    }

    #[test]
    fn maintains_high_precision_on_values() {
        let polar = Degrees::new(90.654_321);   // Dec
        let azimuth = Degrees::new(45.123_456); // RA

        let coord = ICRS::new(polar, azimuth);

        // ra() = azimuth, dec() = polar
        assert!((coord.ra().value() - 45.123_456).abs() < 1e-6);
        assert!((coord.dec().value() - 90.654_321).abs() < 1e-6);
    }

    const EPS: f64 = 1e-6;

    #[test]
    fn position_method_promotes_with_given_radius() {
        use crate::coordinates::centers::Barycentric;

        // Using new_icrs(ra, dec) for correct ICRS convention
        let dir: ICRS = ICRS::new_icrs(Degrees::new(120.0), Degrees::new(-30.0));
        let pos = dir.position::<Barycentric, AstronomicalUnit>(Quantity::<AstronomicalUnit>::new(2.0));

        // angles are preserved (RA=120, Dec=-30)
        assert!((pos.ra().value() - 120.0).abs() < EPS, "RA mismatch: got {}", pos.ra().value());
        assert!((pos.dec().value() - (-30.0)).abs() < EPS, "Dec mismatch: got {}", pos.dec().value());

        // distance matches the supplied magnitude
        assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
    }

    #[test]
    fn direction_display_mentions_frame() {
        let eq: Equatorial = Equatorial::new(Degrees::new(45.0), Degrees::new(10.0));
        let s = eq.to_string();
        assert!(s.contains("Equatorial"), "missing frame");
        // No center in directions anymore
        assert!(!s.contains("Center:"), "should not have center");
    }

    #[test]
    fn angular_separation_identity() {
        let a = ICRS::new(Degrees::new(45.0), Degrees::new(30.0));
        let sep = a.angular_separation(&a);
        assert!(sep.abs().value() < 1e-10, "expected 0°, got {}", sep);
    }
}
