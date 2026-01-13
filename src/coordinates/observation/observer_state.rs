//! # Observer State
//!
//! Encapsulates the complete state of an observer required for computing
//! observation-model effects like aberration.
//!
//! The observer state includes:
//! - Position of the observer (for parallax)
//! - Velocity of the observer (for aberration)
//! - The epoch of observation
//!
//! ## Why Observer State?
//!
//! Aberration depends on the observer's velocity relative to the observed light rays.
//! This velocity has two main components:
//!
//! 1. **Earth's orbital motion** around the Sun (annual aberration, ~20.5")
//! 2. **Observer's diurnal motion** due to Earth's rotation (diurnal aberration, ~0.3")
//!
//! By encapsulating this in `ObserverState`, we ensure that aberration cannot be
//! applied without explicit observer information.

use crate::astro::JulianDate;
use crate::bodies::solar_system::Earth;
use crate::coordinates::cartesian::Velocity;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::frames::EquatorialMeanJ2000;
use qtty::{AstronomicalUnit, Day};

/// Velocity unit: AU per day
pub type AuPerDay = qtty::Per<AstronomicalUnit, Day>;

/// The complete state of an observer required for observation-model effects.
///
/// This struct captures everything needed to compute aberration and other
/// observer-dependent corrections:
///
/// - **velocity**: The observer's velocity (for aberration)
/// - **jd**: The Julian Date of observation
///
/// # Velocity Components
///
/// For a ground-based observer, the velocity includes:
/// - Earth's orbital velocity around the Sun
/// - (Future) Diurnal velocity from Earth's rotation
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::observation::ObserverState;
/// use siderust::astro::JulianDate;
///
/// // Create observer state for a geocentric observer
/// let obs = ObserverState::geocentric(JulianDate::J2000);
///
/// // The observer velocity is now available for aberration calculations
/// println!("Observer velocity: {:?}", obs.velocity());
/// ```
#[derive(Debug, Clone)]
pub struct ObserverState {
    /// Observer velocity in equatorial coordinates (AU/day)
    velocity: Velocity<EquatorialMeanJ2000, AuPerDay>,
    /// Julian Date of observation
    jd: JulianDate,
}

impl ObserverState {
    /// Creates an observer state for an observer at Earth's center.
    ///
    /// This is the standard case for computing annual aberration. The velocity
    /// is Earth's heliocentric orbital velocity.
    ///
    /// # Arguments
    ///
    /// * `jd` - The Julian Date of observation
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::observation::ObserverState;
    /// use siderust::astro::JulianDate;
    ///
    /// let obs = ObserverState::geocentric(JulianDate::J2000);
    /// ```
    pub fn geocentric(jd: JulianDate) -> Self {
        use crate::coordinates::transform::TransformFrame;

        // Get Earth's heliocentric velocity from VSOP87A
        let vel_ecl = Earth::vsop87a_vel(jd);

        // Transform to equatorial frame
        let velocity: Velocity<EquatorialMeanJ2000, AuPerDay> = vel_ecl.to_frame();

        Self { velocity, jd }
    }

    /// Creates an observer state for a topocentric observer.
    ///
    /// This includes Earth's orbital velocity plus the diurnal velocity
    /// from the observer's location on Earth's surface.
    ///
    /// # Arguments
    ///
    /// * `site` - The observer's location on Earth
    /// * `jd` - The Julian Date of observation
    ///
    /// # Note
    ///
    /// Currently this only includes Earth's orbital velocity (annual aberration).
    /// Diurnal aberration (~0.3") is not yet implemented.
    pub fn topocentric(site: &ObserverSite, jd: JulianDate) -> Self {
        use crate::coordinates::transform::TransformFrame;

        // Get Earth's heliocentric velocity from VSOP87A
        let vel_ecl = Earth::vsop87a_vel(jd);

        // Transform to equatorial frame
        let velocity: Velocity<EquatorialMeanJ2000, AuPerDay> = vel_ecl.to_frame();

        // TODO: Add diurnal velocity from Earth rotation
        // For a complete implementation, we would add the observer's
        // velocity due to Earth's rotation (requires GMST and site position)
        let _ = site; // Suppress unused warning for now

        Self { velocity, jd }
    }

    /// Creates an observer state from explicit velocity components.
    ///
    /// Use this for custom observer scenarios (e.g., spacecraft).
    ///
    /// # Arguments
    ///
    /// * `velocity` - Observer velocity in equatorial coordinates
    /// * `jd` - The Julian Date of observation
    pub fn from_velocity(
        velocity: Velocity<EquatorialMeanJ2000, AuPerDay>,
        jd: JulianDate,
    ) -> Self {
        Self { velocity, jd }
    }

    /// Returns the observer's velocity in equatorial coordinates.
    pub fn velocity(&self) -> &Velocity<EquatorialMeanJ2000, AuPerDay> {
        &self.velocity
    }

    /// Returns the Julian Date of observation.
    pub fn jd(&self) -> JulianDate {
        self.jd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geocentric_observer_state() {
        let jd = JulianDate::J2000;
        let obs = ObserverState::geocentric(jd);

        // Earth's orbital velocity is approximately 30 km/s
        // In AU/day: ~30 km/s * 86400 s/day / 149597870.7 km/AU ≈ 0.017 AU/day
        let vel = obs.velocity();
        let speed =
            (vel.x().value().powi(2) + vel.y().value().powi(2) + vel.z().value().powi(2)).sqrt();

        // Speed should be around 0.017 AU/day (Earth's orbital velocity)
        assert!(speed > 0.015, "Earth orbital speed too low: {}", speed);
        assert!(speed < 0.020, "Earth orbital speed too high: {}", speed);
    }

    #[test]
    fn test_observer_jd() {
        let jd = JulianDate::new(2451545.0);
        let obs = ObserverState::geocentric(jd);
        assert_eq!(obs.jd().value(), 2451545.0);
    }
}
