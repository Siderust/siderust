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

use crate::astro::sidereal::unmodded_gst;
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

        // Use SSB-referenced (barycentric) Earth velocity for annual aberration.
        let vel_ecl = Earth::vsop87e_vel(jd);

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
    /// Diurnal aberration (~0.3") is included via an Earth-rotation model based on GMST.
    pub fn topocentric(site: &ObserverSite, jd: JulianDate) -> Self {
        use crate::coordinates::transform::TransformFrame;
        use qtty::{Meter, Radian};

        // Annual (orbital) component: barycentric Earth velocity (VSOP87E).
        let vel_ecl = Earth::vsop87e_vel(jd);

        // Transform to equatorial frame
        let mut velocity: Velocity<EquatorialMeanJ2000, AuPerDay> = vel_ecl.to_frame();

        // Diurnal (rotational) component: v = ω × r, computed in ECEF then rotated to equatorial.
        // This uses GMST as a first-order Earth rotation model (UT1 should be supplied via `jd`).
        const AU_M: f64 = 149_597_870_700.0;
        const SECONDS_PER_DAY: f64 = 86_400.0;
        // IERS Conventions 2010/2020: Earth rotation rate (rad/s), nominal.
        const OMEGA_EARTH: f64 = 7.292_115_0e-5;

        let site_itrf_m = site.geocentric_itrf::<Meter>();
        let rx = site_itrf_m.x().value();
        let ry = site_itrf_m.y().value();

        // ω = (0,0,OMEGA_EARTH) in ECEF => ω×r = (-ω*y, ω*x, 0)
        let vx_ecef_mps = -OMEGA_EARTH * ry;
        let vy_ecef_mps = OMEGA_EARTH * rx;
        let vz_ecef_mps = 0.0;

        // Rotate ECEF velocity into the mean equator/equinox of J2000 using GMST about +Z.
        let gmst_rad = unmodded_gst(jd).to::<Radian>().value();
        let (sin_g, cos_g) = gmst_rad.sin_cos();

        let vx_eq_mps = vx_ecef_mps * cos_g - vy_ecef_mps * sin_g;
        let vy_eq_mps = vx_ecef_mps * sin_g + vy_ecef_mps * cos_g;
        let vz_eq_mps = vz_ecef_mps;

        let mps_to_au_per_day = |v_mps: f64| v_mps * SECONDS_PER_DAY / AU_M;
        let v_diurnal = Velocity::<EquatorialMeanJ2000, AuPerDay>::new(
            qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(mps_to_au_per_day(vx_eq_mps)),
            qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(mps_to_au_per_day(vy_eq_mps)),
            qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(mps_to_au_per_day(vz_eq_mps)),
        );

        velocity = Velocity::<EquatorialMeanJ2000, AuPerDay>::new(
            velocity.x() + v_diurnal.x(),
            velocity.y() + v_diurnal.y(),
            velocity.z() + v_diurnal.z(),
        );

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
    use qtty::*;

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

    #[test]
    fn test_topocentric_includes_diurnal_velocity() {
        let jd = JulianDate::J2000;
        let geo = ObserverState::geocentric(jd);
        let site = ObserverSite::new(0.0 * DEG, 0.0 * DEG, 0.0 * M); // equator
        let topo = ObserverState::topocentric(&site, jd);

        let dvx = topo.velocity().x().value() - geo.velocity().x().value();
        let dvy = topo.velocity().y().value() - geo.velocity().y().value();
        let dvz = topo.velocity().z().value() - geo.velocity().z().value();
        let dv = (dvx * dvx + dvy * dvy + dvz * dvz).sqrt();

        // Equatorial surface speed is ~465 m/s ≈ 2.685e-4 AU/day.
        assert!(dv > 2.3e-4, "diurnal speed too low: {}", dv);
        assert!(dv < 3.1e-4, "diurnal speed too high: {}", dv);
    }
}
