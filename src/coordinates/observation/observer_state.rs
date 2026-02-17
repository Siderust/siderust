// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

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

use crate::astro::eop::EopProvider;
use crate::calculus::ephemeris::Ephemeris;
use crate::coordinates::cartesian::Velocity;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::frames::EquatorialMeanJ2000;
use crate::coordinates::transform::centers::position::to_topocentric::itrs_to_equatorial_mean_j2000_rotation;
use crate::coordinates::transform::context::{AstroContext, DefaultEphemeris};
use crate::time::JulianDate;
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
/// use siderust::time::JulianDate;
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
    /// use siderust::time::JulianDate;
    ///
    /// let obs = ObserverState::geocentric(JulianDate::J2000);
    /// ```
    pub fn geocentric(jd: JulianDate) -> Self {
        use crate::coordinates::transform::TransformFrame;

        // Use SSB-referenced (barycentric) Earth velocity for annual aberration.
        let vel_ecl = DefaultEphemeris::earth_barycentric_velocity(jd);

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
    /// Diurnal aberration (~0.3") uses the default EOP-backed Earth-rotation chain.
    pub fn topocentric(site: &ObserverSite, jd: JulianDate) -> Self {
        let ctx: AstroContext = AstroContext::default();
        Self::topocentric_with_ctx(site, jd, &ctx)
    }

    /// Context-aware topocentric observer state.
    pub fn topocentric_with_ctx<Eph, Eop: EopProvider, Nut>(
        site: &ObserverSite,
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Self {
        use crate::coordinates::transform::TransformFrame;
        use qtty::{Meter, Second, Seconds};

        // Annual (orbital) component: barycentric Earth velocity.
        let vel_ecl = DefaultEphemeris::earth_barycentric_velocity(jd);

        // Transform to equatorial frame
        let mut velocity: Velocity<EquatorialMeanJ2000, AuPerDay> = vel_ecl.to_frame();

        // Diurnal (rotational) component: v = ω × r, computed in ECEF then
        // rotated through the same high-precision terrestrial->celestial chain
        // used by topocentric position transforms.
        // IERS Conventions 2010/2020: Earth rotation rate (rad/s), nominal.
        const OMEGA_EARTH: f64 = 7.292_115_0e-5;

        let site_itrf_m = site.geocentric_itrf::<Meter>();
        type MetersPerSecond = qtty::Per<Meter, Second>;
        let one_second = Seconds::new(1.0);

        // ω = (0,0,OMEGA_EARTH) in ECEF => ω×r = (-ω*y, ω*x, 0)
        let vx_ecef: qtty::Quantity<MetersPerSecond> =
            (-site_itrf_m.y() * OMEGA_EARTH) / one_second;
        let vy_ecef: qtty::Quantity<MetersPerSecond> = (site_itrf_m.x() * OMEGA_EARTH) / one_second;
        let vz_ecef: qtty::Quantity<MetersPerSecond> = qtty::Quantity::zero();

        let rot = itrs_to_equatorial_mean_j2000_rotation(jd, ctx);
        let [vx_eq, vy_eq, vz_eq] = rot * [vx_ecef, vy_ecef, vz_ecef];

        let v_diurnal = Velocity::<EquatorialMeanJ2000, AuPerDay>::new(
            vx_eq.to::<AuPerDay>(),
            vy_eq.to::<AuPerDay>(),
            vz_eq.to::<AuPerDay>(),
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
        let speed = (vel.x() * vel.x() + vel.y() * vel.y() + vel.z() * vel.z()).sqrt();

        // Speed should be around 0.017 AU/day (Earth's orbital velocity)
        assert!(speed > 0.015, "Earth orbital speed too low: {}", speed);
        assert!(speed < 0.020, "Earth orbital speed too high: {}", speed);
    }

    #[test]
    fn test_observer_jd() {
        let jd = JulianDate::new(2451545.0);
        let obs = ObserverState::geocentric(jd);
        assert_eq!(obs.jd(), jd);
    }

    #[test]
    fn test_topocentric_includes_diurnal_velocity() {
        let jd = JulianDate::J2000;
        let geo = ObserverState::geocentric(jd);
        let site = ObserverSite::new(0.0 * DEG, 0.0 * DEG, 0.0 * M); // equator
        let topo = ObserverState::topocentric(&site, jd);

        let dvx = topo.velocity().x() - geo.velocity().x();
        let dvy = topo.velocity().y() - geo.velocity().y();
        let dvz = topo.velocity().z() - geo.velocity().z();
        let dv = (dvx * dvx + dvy * dvy + dvz * dvz).sqrt();

        // Equatorial surface speed is ~465 m/s ≈ 2.685e-4 AU/day.
        assert!(dv > 2.3e-4, "diurnal speed too low: {}", dv);
        assert!(dv < 3.1e-4, "diurnal speed too high: {}", dv);
    }
}
