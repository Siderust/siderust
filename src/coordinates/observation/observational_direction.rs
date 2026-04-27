// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observational Direction Types
//!
//! This module defines direction types with explicit observational state markers.
//! A direction can be either `Astrometric` (geometric) or `Apparent` (observed).
//!
//! ## Type-Level Observational State
//!
//! Rather than storing state as runtime metadata, we use wrapper types:
//!
//! - `Astrometric<D>`: A geometric direction, not corrected for observer motion
//! - `Apparent<D>`: A direction as observed, with aberration applied
//!
//! This makes the physical meaning unambiguous at compile time.
//!
//! ## Conversions
//!
//! Converting between states requires an `ObserverState`:
//!
//! ```rust
//! use siderust::coordinates::observation::{Astrometric, Apparent, ObserverState};
//! use siderust::coordinates::spherical::direction::EquatorialMeanJ2000;
//! use siderust::time::JulianDate;
//! use siderust::qtty::*;
//!
//! let geo_dir = Astrometric::new(EquatorialMeanJ2000::new(45.0 * DEG, 20.0 * DEG));
//! let obs = ObserverState::geocentric(JulianDate::J2000);
//!
//! // Explicit conversion with observer state
//! let app_dir: Apparent<EquatorialMeanJ2000> = geo_dir.to_apparent(&obs);
//!
//! // Can convert back
//! let geo_again: Astrometric<EquatorialMeanJ2000> = app_dir.to_astrometric(&obs);
//! ```

use super::ObserverState;
use crate::astro::aberration::{
    apply_aberration_to_direction_with_velocity, remove_aberration_from_direction_with_velocity,
};
use crate::coordinates::cartesian;
use crate::coordinates::frames::{EquatorialMeanJ2000, MutableFrame};
use crate::coordinates::spherical;
use crate::coordinates::transform::TransformFrame;

// =============================================================================
// Astrometric Direction
// =============================================================================

/// A direction representing the geometric (astrometric) position of a celestial object.
///
/// Astrometric directions are:
/// - Corrected for proper motion and parallax (if applicable)
/// - **Not** corrected for aberration or other observer-motion effects
/// - The "true" direction to the object in space
///
/// To obtain an apparent direction (as observed), call `to_apparent()` with
/// an `ObserverState`.
///
/// # Type Parameter
///
/// - `D`: The underlying direction type (e.g., `spherical::direction::EquatorialMeanJ2000`)
#[derive(Debug, Clone, Copy)]
pub struct Astrometric<D> {
    direction: D,
}

impl<D> Astrometric<D> {
    /// Creates a new astrometric direction from a geometric direction.
    ///
    /// # Arguments
    ///
    /// * `direction` - The geometric direction
    pub fn new(direction: D) -> Self {
        Self { direction }
    }

    /// Returns a reference to the underlying direction.
    pub fn direction(&self) -> &D {
        &self.direction
    }

    /// Consumes the wrapper and returns the underlying direction.
    pub fn into_inner(self) -> D {
        self.direction
    }
}

impl<F: MutableFrame> Astrometric<spherical::Direction<F>> {
    /// Converts this astrometric direction to an apparent direction.
    ///
    /// This applies aberration based on the observer's velocity.
    ///
    /// # Arguments
    ///
    /// * `obs` - The observer state (provides velocity for aberration)
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::observation::{Astrometric, Apparent, ObserverState};
    /// use siderust::coordinates::spherical::direction::EquatorialMeanJ2000;
    /// use siderust::time::JulianDate;
    /// use siderust::qtty::*;
    ///
    /// let astrometric = Astrometric::new(EquatorialMeanJ2000::new(0.0 * DEG, 0.0 * DEG));
    /// let obs = ObserverState::geocentric(JulianDate::J2000);
    /// let apparent: Apparent<EquatorialMeanJ2000> = astrometric.to_apparent(&obs);
    /// ```
    pub fn to_apparent(self, obs: &ObserverState) -> Apparent<spherical::Direction<F>>
    where
        cartesian::Direction<EquatorialMeanJ2000>: TransformFrame<cartesian::Direction<F>>,
        cartesian::Direction<F>: TransformFrame<cartesian::Direction<EquatorialMeanJ2000>>,
    {
        // Convert to cartesian for vector operations
        let cart_dir = self.direction.to_cartesian();

        // Transform to equatorial frame for aberration
        let cart_eq: cartesian::Direction<EquatorialMeanJ2000> = cart_dir.to_frame();

        // Get observer velocity
        let vel = obs.velocity();

        let apparent_eq = apply_aberration_to_direction_with_velocity(cart_eq, vel);

        // Transform back to original frame
        let apparent_cart: cartesian::Direction<F> = apparent_eq.to_frame();
        let apparent_sph = apparent_cart.to_spherical();

        Apparent::new(apparent_sph)
    }
}

impl<F: MutableFrame> Astrometric<cartesian::Direction<F>> {
    /// Converts this astrometric direction to an apparent direction.
    ///
    /// This applies aberration based on the observer's velocity.
    pub fn to_apparent(self, obs: &ObserverState) -> Apparent<cartesian::Direction<F>>
    where
        cartesian::Direction<EquatorialMeanJ2000>: TransformFrame<cartesian::Direction<F>>,
        cartesian::Direction<F>: TransformFrame<cartesian::Direction<EquatorialMeanJ2000>>,
    {
        // Transform to equatorial frame for aberration
        let cart_eq: cartesian::Direction<EquatorialMeanJ2000> = self.direction.to_frame();

        // Get observer velocity
        let vel = obs.velocity();

        let apparent_eq = apply_aberration_to_direction_with_velocity(cart_eq, vel);

        // Transform back to original frame
        let apparent_cart: cartesian::Direction<F> = apparent_eq.to_frame();

        Apparent::new(apparent_cart)
    }
}

// =============================================================================
// Apparent Direction
// =============================================================================

/// A direction representing the apparent (observed) position of a celestial object.
///
/// Apparent directions include:
/// - All astrometric corrections (proper motion, parallax)
/// - Aberration correction based on observer velocity
/// - (Future) Other observation effects like gravitational deflection
///
/// The direction has been shifted from its geometric (astrometric) position
/// due to the finite speed of light and the observer's motion.
///
/// # Type Parameter
///
/// - `D`: The underlying direction type (e.g., `spherical::direction::EquatorialMeanJ2000`)
#[derive(Debug, Clone, Copy)]
pub struct Apparent<D> {
    direction: D,
}

impl<D> Apparent<D> {
    /// Creates a new apparent direction.
    ///
    /// # Note
    ///
    /// You typically obtain apparent directions by calling `to_apparent()`
    /// on an `Astrometric` direction. This constructor is provided for
    /// cases where you have a known apparent position.
    pub fn new(direction: D) -> Self {
        Self { direction }
    }

    /// Returns a reference to the underlying direction.
    pub fn direction(&self) -> &D {
        &self.direction
    }

    /// Consumes the wrapper and returns the underlying direction.
    pub fn into_inner(self) -> D {
        self.direction
    }
}

impl<F: MutableFrame> Apparent<spherical::Direction<F>> {
    /// Converts this apparent direction back to an astrometric direction.
    ///
    /// This removes aberration based on the observer's velocity.
    ///
    /// # Arguments
    ///
    /// * `obs` - The observer state (provides velocity for aberration removal)
    pub fn to_astrometric(self, obs: &ObserverState) -> Astrometric<spherical::Direction<F>>
    where
        cartesian::Direction<EquatorialMeanJ2000>: TransformFrame<cartesian::Direction<F>>,
        cartesian::Direction<F>: TransformFrame<cartesian::Direction<EquatorialMeanJ2000>>,
    {
        // Convert to cartesian for vector operations
        let cart_dir = self.direction.to_cartesian();

        // Transform to equatorial frame for aberration
        let cart_eq: cartesian::Direction<EquatorialMeanJ2000> = cart_dir.to_frame();

        // Get observer velocity
        let vel = obs.velocity();

        let astrometric_eq = remove_aberration_from_direction_with_velocity(cart_eq, vel);

        // Transform back to original frame
        let astrometric_cart: cartesian::Direction<F> = astrometric_eq.to_frame();
        let astrometric_sph = astrometric_cart.to_spherical();

        Astrometric::new(astrometric_sph)
    }
}

impl<F: MutableFrame> Apparent<cartesian::Direction<F>> {
    /// Converts this apparent direction back to an astrometric direction.
    ///
    /// This removes aberration based on the observer's velocity.
    pub fn to_astrometric(self, obs: &ObserverState) -> Astrometric<cartesian::Direction<F>>
    where
        cartesian::Direction<EquatorialMeanJ2000>: TransformFrame<cartesian::Direction<F>>,
        cartesian::Direction<F>: TransformFrame<cartesian::Direction<EquatorialMeanJ2000>>,
    {
        // Transform to equatorial frame for aberration
        let cart_eq: cartesian::Direction<EquatorialMeanJ2000> = self.direction.to_frame();

        // Get observer velocity
        let vel = obs.velocity();

        let astrometric_eq = remove_aberration_from_direction_with_velocity(cart_eq, vel);

        // Transform back to original frame
        let astrometric_cart: cartesian::Direction<F> = astrometric_eq.to_frame();

        Astrometric::new(astrometric_cart)
    }
}

// =============================================================================
// Display implementations
// =============================================================================

impl<D: std::fmt::Display> std::fmt::Display for Astrometric<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Astrometric({})", self.direction)
    }
}

impl<D: std::fmt::Display> std::fmt::Display for Apparent<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Apparent({})", self.direction)
    }
}

impl<D: std::fmt::LowerExp> std::fmt::LowerExp for Astrometric<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Astrometric(")?;
        std::fmt::LowerExp::fmt(&self.direction, f)?;
        write!(f, ")")
    }
}

impl<D: std::fmt::LowerExp> std::fmt::LowerExp for Apparent<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Apparent(")?;
        std::fmt::LowerExp::fmt(&self.direction, f)?;
        write!(f, ")")
    }
}

impl<D: std::fmt::UpperExp> std::fmt::UpperExp for Astrometric<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Astrometric(")?;
        std::fmt::UpperExp::fmt(&self.direction, f)?;
        write!(f, ")")
    }
}

impl<D: std::fmt::UpperExp> std::fmt::UpperExp for Apparent<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Apparent(")?;
        std::fmt::UpperExp::fmt(&self.direction, f)?;
        write!(f, ")")
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::JulianDate;
    use crate::qtty::*;

    #[test]
    fn test_astrometric_to_apparent_introduces_shift() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        // Create an astrometric direction
        let astrometric = Astrometric::new(spherical::direction::EquatorialMeanJ2000::new(
            0.0 * DEG,
            0.0 * DEG,
        ));

        // Convert to apparent
        let apparent: Apparent<spherical::direction::EquatorialMeanJ2000> =
            astrometric.to_apparent(&obs);

        // The shift should be on the order of 20 arcseconds (aberration constant)
        let original = astrometric.direction();
        let shifted = apparent.direction();

        // Handle azimuth wrap-around (e.g., 359.999° vs 0.001° is actually 0.002° apart)
        let mut delta_ra = (shifted.azimuth - original.azimuth).abs();
        if delta_ra > Degrees::new(180.0) {
            delta_ra = Degrees::new(360.0) - delta_ra;
        }
        let delta_dec = (shifted.polar - original.polar).abs();

        // At least one should have changed
        assert!(
            delta_ra > Degrees::new(0.0) || delta_dec > Degrees::new(0.0),
            "Expected aberration to introduce a shift"
        );

        // Shift should be small (less than 1 degree, typically ~20 arcsec = 0.006 deg)
        assert!(
            delta_ra < Degrees::new(0.1) && delta_dec < Degrees::new(0.1),
            "Aberration shift too large: dRA={}, dDec={}",
            delta_ra,
            delta_dec
        );
    }

    #[test]
    fn test_roundtrip_preserves_direction() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        // Create an astrometric direction
        let original = Astrometric::new(spherical::direction::EquatorialMeanJ2000::new(
            45.0 * DEG,
            30.0 * DEG,
        ));

        // Convert to apparent and back
        let apparent = original.to_apparent(&obs);
        let recovered: Astrometric<spherical::direction::EquatorialMeanJ2000> =
            apparent.to_astrometric(&obs);

        // Should be very close to original
        let orig_dir = original.direction();
        let rec_dir = recovered.direction();

        let delta_ra = (rec_dir.azimuth - orig_dir.azimuth).abs();
        let delta_dec = (rec_dir.polar - orig_dir.polar).abs();

        // Tolerance of 1e-6 degrees ≈ 0.003 arcseconds (numerical precision limit)
        assert!(
            delta_ra < Degrees::new(1e-6),
            "RA not preserved in roundtrip: delta = {}",
            delta_ra
        );
        assert!(
            delta_dec < Degrees::new(1e-6),
            "Dec not preserved in roundtrip: delta = {}",
            delta_dec
        );
    }

    #[test]
    fn test_apparent_at_pole() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        // At the north pole
        let astrometric = Astrometric::new(spherical::direction::EquatorialMeanJ2000::new(
            0.0 * DEG,
            90.0 * DEG,
        ));

        let apparent = astrometric.to_apparent(&obs);

        // Declination should decrease slightly
        assert!(
            apparent.direction().polar < Degrees::new(90.0),
            "Declination at pole should decrease due to aberration"
        );
    }

    #[test]
    fn test_type_safety() {
        // This test verifies that you can't accidentally mix astrometric and apparent
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        let astrometric = Astrometric::new(spherical::direction::EquatorialMeanJ2000::new(
            10.0 * DEG,
            20.0 * DEG,
        ));
        let apparent = astrometric.to_apparent(&obs);

        // These are different types - can't be confused
        let _astrometric_inner: &spherical::direction::EquatorialMeanJ2000 =
            astrometric.direction();
        let _apparent_inner: &spherical::direction::EquatorialMeanJ2000 = apparent.direction();

        // The wrappers make the distinction explicit
        assert!(format!("{}", astrometric).contains("Astrometric"));
        assert!(format!("{}", apparent).contains("Apparent"));
    }

    // =====================================================================
    // Cartesian direction paths
    // =====================================================================

    #[test]
    fn test_cartesian_astrometric_to_apparent_roundtrip() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        let cart_dir = cartesian::Direction::<crate::coordinates::frames::EquatorialMeanJ2000>::new(
            0.6, 0.7, 0.3,
        );
        let astrometric = Astrometric::new(cart_dir);

        let apparent = astrometric.to_apparent(&obs);
        let recovered = apparent.to_astrometric(&obs);

        let orig = astrometric.direction();
        let rec = recovered.direction();
        assert!((rec.x() - orig.x()).abs() < 1e-8);
        assert!((rec.y() - orig.y()).abs() < 1e-8);
        assert!((rec.z() - orig.z()).abs() < 1e-8);
    }

    #[test]
    fn test_cartesian_aberration_introduces_shift() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        let cart_dir = cartesian::Direction::<crate::coordinates::frames::EquatorialMeanJ2000>::new(
            1.0, 0.0, 0.0,
        );
        let astrometric = Astrometric::new(cart_dir);
        let apparent = astrometric.to_apparent(&obs);

        let orig = astrometric.direction();
        let shifted = apparent.direction();
        let delta = (shifted.x() - orig.x()).abs()
            + (shifted.y() - orig.y()).abs()
            + (shifted.z() - orig.z()).abs();
        assert!(delta > 0.0, "Aberration should shift the direction");
        assert!(
            delta < 0.001,
            "Aberration shift should be small, got {}",
            delta
        );
    }

    #[test]
    fn test_cartesian_into_inner() {
        let cart_dir = cartesian::Direction::<crate::coordinates::frames::EquatorialMeanJ2000>::new(
            0.0, 0.0, 1.0,
        );
        let astrometric = Astrometric::new(cart_dir);
        let inner = astrometric.into_inner();
        assert!((inner.z() - 1.0).abs() < 1e-12);
    }

    // =====================================================================
    // Display / LowerExp / UpperExp formatting
    // =====================================================================

    #[test]
    fn test_display_formatting() {
        let dir = spherical::direction::EquatorialMeanJ2000::new(10.0 * DEG, 20.0 * DEG);
        let astrometric = Astrometric::new(dir);
        let display = format!("{}", astrometric);
        assert!(display.starts_with("Astrometric("));
        assert!(display.ends_with(')'));
    }

    #[test]
    fn test_apparent_display_formatting() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);
        let astrometric = Astrometric::new(spherical::direction::EquatorialMeanJ2000::new(
            10.0 * DEG,
            20.0 * DEG,
        ));
        let apparent = astrometric.to_apparent(&obs);
        let display = format!("{}", apparent);
        assert!(display.starts_with("Apparent("));
        assert!(display.ends_with(')'));
    }

    #[test]
    fn test_lower_exp_formatting() {
        let dir = spherical::direction::EquatorialMeanJ2000::new(10.0 * DEG, 20.0 * DEG);
        let astrometric = Astrometric::new(dir);
        let display = format!("{:e}", astrometric);
        assert!(display.starts_with("Astrometric("));

        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);
        let apparent = astrometric.to_apparent(&obs);
        let display = format!("{:e}", apparent);
        assert!(display.starts_with("Apparent("));
    }

    #[test]
    fn test_upper_exp_formatting() {
        let dir = spherical::direction::EquatorialMeanJ2000::new(10.0 * DEG, 20.0 * DEG);
        let astrometric = Astrometric::new(dir);
        let display = format!("{:E}", astrometric);
        assert!(display.starts_with("Astrometric("));

        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);
        let apparent = astrometric.to_apparent(&obs);
        let display = format!("{:E}", apparent);
        assert!(display.starts_with("Apparent("));
    }
}
