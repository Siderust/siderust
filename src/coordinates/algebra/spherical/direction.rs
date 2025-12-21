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
//! ## Storage
//!
//! A spherical direction stores only the two angles (polar, azimuth). The implicit
//! radius is always 1. This contrasts with [`Position`](super::Position), which
//! stores a radial distance in addition to the angles.
//!
//! ## Observer-Dependent Directions
//!
//! To compute the direction from an observer to a target (line of sight), use
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) with two
//! positions. This correctly accounts for the observer's location.
//!
//! ## Angle Convention
//!
//! - **polar (θ)**: Latitude / declination / altitude — range `[-90°, +90°]`
//! - **azimuth (φ)**: Longitude / right ascension / azimuth — normalized to `[0°, 360°)`
//!
//! Constructors canonicalize angles to these ranges.
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
//! // Direction to Vega (J2000) using frame-specific constructor
//! // Note: new_icrs takes (RA, Dec) in astronomical convention
//! let vega: ICRS = ICRS::new_icrs(279.23473479 * DEG, 38.78368896 * DEG);
//!
//! // Convert that direction into a position one parsec away:
//! use siderust::coordinates::algebra::centers::Barycentric;
//! let one_pc = 1.0 * PC;
//! let pos = vega.position::<Barycentric, _>(one_pc);
//! println!("Vega direction = {vega}\nVega@1pc position = {pos}");
//! ```

use crate::coordinates::algebra::centers::ReferenceCenter;
use crate::coordinates::algebra::frames::{self, ReferenceFrame};
use qtty::{Degrees, LengthUnit, Quantity};

use std::marker::PhantomData;

// =============================================================================
// Angle Canonicalization Helpers
// =============================================================================

/// Normalizes an azimuth angle to the range `[0°, 360°)`.
#[inline]
fn canonicalize_azimuth(azimuth: Degrees) -> Degrees {
    let mut val = azimuth.value() % 360.0;
    if val < 0.0 {
        val += 360.0;
    }
    Degrees::new(val)
}

/// Clamps a polar angle to `[-90°, +90°]`.
///
/// This is the mathematically correct range for latitude/declination.
/// Values outside this range are clamped (not wrapped).
#[inline]
fn canonicalize_polar(polar: Degrees) -> Degrees {
    Degrees::new(polar.value().clamp(-90.0, 90.0))
}

/// A spherical direction (unit vector) in a specific reference frame.
///
/// Directions are frame-dependent but center-independent (free vectors).
/// They cannot undergo center transformations, only frame transformations.
///
/// The direction stores only angles; the implicit radius is always 1.
/// For a spherical coordinate with explicit distance, use [`Position`](super::Position).
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
///
/// # Invariants
///
/// - `polar` is in `[-90°, +90°]`
/// - `azimuth` is in `[0°, 360°)`
#[derive(Debug, Clone, Copy)]
pub struct Direction<F: ReferenceFrame> {
    /// Polar angle (θ) - latitude, declination, or altitude, in degrees.
    /// Range: `[-90°, +90°]`
    pub polar: Degrees,
    /// Azimuthal angle (φ) - longitude, right ascension, or azimuth, in degrees.
    /// Range: `[0°, 360°)`
    pub azimuth: Degrees,
    _frame: PhantomData<F>,
}

impl<F: ReferenceFrame> Direction<F> {
    /// Creates a new direction from polar and azimuth angles.
    ///
    /// Angles are canonicalized:
    /// - `polar` is clamped to `[-90°, +90°]`
    /// - `azimuth` is normalized to `[0°, 360°)`
    ///
    /// * `polar`: angle from the reference plane (latitude, declination), in degrees
    /// * `azimuth`: angle from the reference direction (longitude, RA), in degrees
    ///
    /// # Note
    ///
    /// For frame-specific constructors with more intuitive parameter order
    /// (e.g., `new_icrs(ra, dec)`, `new_ecliptic(lon, lat)`), see the extension
    /// methods in [`crate::coordinates::astro::spherical`].
    pub fn new(polar: Degrees, azimuth: Degrees) -> Self {
        Self {
            polar: canonicalize_polar(polar),
            azimuth: canonicalize_azimuth(azimuth),
            _frame: PhantomData,
        }
    }

    /// Creates a direction from raw angle values without canonicalization.
    ///
    /// # Safety
    ///
    /// The caller must ensure angles are within canonical ranges:
    /// - `polar` in `[-90°, +90°]`
    /// - `azimuth` in `[0°, 360°)`
    ///
    /// Use [`new`](Self::new) for automatic canonicalization.
    pub const fn new_raw(polar: Degrees, azimuth: Degrees) -> Self {
        Self {
            polar,
            azimuth,
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
        use qtty::{Degree, Radian, Radians};

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
    pub fn to_cartesian(&self) -> crate::coordinates::algebra::cartesian::Direction<F>
    where
        F: frames::MutableFrame,
    {
        use qtty::Radian;

        let polar_rad = self.polar.to::<Radian>();
        let azimuth_rad = self.azimuth.to::<Radian>();

        let x = azimuth_rad.cos() * polar_rad.cos();
        let y = azimuth_rad.sin() * polar_rad.cos();
        let z = polar_rad.sin();

        crate::coordinates::algebra::cartesian::Direction::<F>::new(x, y, z)
    }

    /// Constructs a spherical direction from a cartesian Direction.
    ///
    /// The resulting angles are canonicalized:
    /// - `polar` in `[-90°, +90°]`
    /// - `azimuth` in `[0°, 360°)`
    pub fn from_cartesian(cart: &crate::coordinates::algebra::cartesian::Direction<F>) -> Self
    where
        F: frames::MutableFrame,
    {
        let x = cart.x();
        let y = cart.y();
        let z = cart.z();

        // Clamp z to prevent asin domain errors from floating-point imprecision
        let z_clamped = z.clamp(-1.0, 1.0);
        let polar = Degrees::new(z_clamped.asin().to_degrees());
        let azimuth = Degrees::new(y.atan2(x).to_degrees());

        Self::new(polar, azimuth)
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

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn creates_valid_spherical_direction() {
        let polar = Degrees::new(45.0);
        let azimuth = Degrees::new(90.0);

        let coord = Direction::<frames::ICRS>::new(polar, azimuth);

        assert_eq!(coord.azimuth.value(), 90.0);
        assert_eq!(coord.polar.value(), 45.0);
        // Direction no longer has a distance field - it's implicitly 1
    }

    #[test]
    fn displays_coordinate_as_string_correctly() {
        let coord = Direction::<frames::ICRS>::new(Degrees::new(30.0), Degrees::new(60.0));
        let output = coord.to_string();
        assert!(output.contains("θ: 30"), "Missing polar angle");
        assert!(output.contains("φ: 60"), "Missing azimuth");
    }

    #[test]
    fn maintains_high_precision_on_values() {
        // Note: polar will be clamped to [-90, 90]
        let polar = Degrees::new(45.123_456);
        let azimuth = Degrees::new(90.654_321);

        let coord = Direction::<frames::ICRS>::new(polar, azimuth);

        assert!((coord.polar.value() - 45.123_456).abs() < 1e-6);
        assert!((coord.azimuth.value() - 90.654_321).abs() < 1e-6);
    }

    const EPS: f64 = 1e-6;

    #[test]
    fn position_method_promotes_with_given_radius() {
        use crate::coordinates::algebra::centers::Barycentric;

        let dir = Direction::<frames::ICRS>::new(Degrees::new(-30.0), Degrees::new(120.0));
        let pos =
            dir.position::<Barycentric, AstronomicalUnit>(Quantity::<AstronomicalUnit>::new(2.0));

        // angles are preserved
        assert!(
            (pos.azimuth.value() - 120.0).abs() < EPS,
            "azimuth mismatch: got {}",
            pos.azimuth.value()
        );
        assert!(
            (pos.polar.value() - (-30.0)).abs() < EPS,
            "polar mismatch: got {}",
            pos.polar.value()
        );

        // distance matches the supplied magnitude
        assert!((pos.distance - 2.0 * AU).abs() < EPS * AU);
    }

    #[test]
    fn direction_display_mentions_frame() {
        let eq = Direction::<frames::Equatorial>::new(Degrees::new(45.0), Degrees::new(10.0));
        let s = eq.to_string();
        assert!(s.contains("Equatorial"), "missing frame");
        // No center in directions anymore
        assert!(!s.contains("Center:"), "should not have center");
    }

    #[test]
    fn angular_separation_identity() {
        let a = Direction::<frames::ICRS>::new(Degrees::new(45.0), Degrees::new(30.0));
        let sep = a.angular_separation(&a);
        assert!(sep.abs().value() < 1e-10, "expected 0°, got {}", sep);
    }

    // =============================================================================
    // Canonicalization Tests
    // =============================================================================

    #[test]
    fn canonicalizes_azimuth_to_positive_range() {
        // Negative azimuth
        let dir = Direction::<frames::ICRS>::new(Degrees::new(0.0), Degrees::new(-90.0));
        assert!((dir.azimuth.value() - 270.0).abs() < EPS);

        // Azimuth > 360
        let dir2 = Direction::<frames::ICRS>::new(Degrees::new(0.0), Degrees::new(450.0));
        assert!((dir2.azimuth.value() - 90.0).abs() < EPS);

        // Large negative
        let dir3 = Direction::<frames::ICRS>::new(Degrees::new(0.0), Degrees::new(-720.0));
        assert!(dir3.azimuth.value().abs() < EPS);
    }

    #[test]
    fn clamps_polar_to_valid_range() {
        // Polar > 90 gets clamped
        let dir = Direction::<frames::ICRS>::new(Degrees::new(100.0), Degrees::new(0.0));
        assert!((dir.polar.value() - 90.0).abs() < EPS);

        // Polar < -90 gets clamped
        let dir2 = Direction::<frames::ICRS>::new(Degrees::new(-100.0), Degrees::new(0.0));
        assert!((dir2.polar.value() - (-90.0)).abs() < EPS);
    }

    // =============================================================================
    // Roundtrip Tests
    // =============================================================================

    #[test]
    fn roundtrip_spherical_cartesian_direction() {
        let original = Direction::<frames::ICRS>::new(Degrees::new(45.0), Degrees::new(30.0));
        let cartesian = original.to_cartesian();
        let recovered = Direction::from_cartesian(&cartesian);

        assert!(
            (recovered.polar.value() - original.polar.value()).abs() < EPS,
            "polar mismatch: {} vs {}",
            recovered.polar.value(),
            original.polar.value()
        );
        assert!(
            (recovered.azimuth.value() - original.azimuth.value()).abs() < EPS,
            "azimuth mismatch: {} vs {}",
            recovered.azimuth.value(),
            original.azimuth.value()
        );
    }

    #[test]
    fn roundtrip_at_poles() {
        // North pole
        let north = Direction::<frames::ICRS>::new(Degrees::new(90.0), Degrees::new(0.0));
        let cart_n = north.to_cartesian();
        let recovered_n = Direction::from_cartesian(&cart_n);
        assert!((recovered_n.polar.value() - 90.0).abs() < EPS);

        // South pole
        let south = Direction::<frames::ICRS>::new(Degrees::new(-90.0), Degrees::new(0.0));
        let cart_s = south.to_cartesian();
        let recovered_s = Direction::from_cartesian(&cart_s);
        assert!((recovered_s.polar.value() - (-90.0)).abs() < EPS);
    }

    #[test]
    fn roundtrip_at_azimuth_boundaries() {
        // At azimuth = 0
        let dir0 = Direction::<frames::ICRS>::new(Degrees::new(30.0), Degrees::new(0.0));
        let cart0 = dir0.to_cartesian();
        let rec0 = Direction::from_cartesian(&cart0);
        assert!((rec0.azimuth.value() - 0.0).abs() < EPS);

        // Near azimuth = 360 (should wrap to near 0)
        let dir360 = Direction::<frames::ICRS>::new(Degrees::new(30.0), Degrees::new(359.9));
        let cart360 = dir360.to_cartesian();
        let rec360 = Direction::from_cartesian(&cart360);
        assert!((rec360.azimuth.value() - 359.9).abs() < EPS);
    }
}
