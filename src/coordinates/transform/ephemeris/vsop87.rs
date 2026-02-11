// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # VSOP87 Ephemeris Backend
//!
//! Implements [`BodyEphemeris`] and [`VelocityEphemeris`] using the built-in
//! VSOP87 planetary theory.

use super::{BodyEphemeris, BodyId, VelocityEphemeris};
use crate::astro::JulianDate;
use crate::bodies::solar_system::{Earth, Jupiter, Mars, Mercury, Neptune, Saturn, Sun, Uranus, Venus};

/// The default ephemeris implementation using VSOP87 planetary theory.
///
/// VSOP87 (Variations Séculaires des Orbites Planétaires) provides high-precision
/// planetary positions for the major planets. This implementation uses:
///
/// - **VSOP87E**: Barycentric rectangular coordinates (ecliptic J2000.0)
/// - **VSOP87A**: Heliocentric rectangular coordinates (ecliptic J2000.0)
///
/// ## Accuracy
///
/// VSOP87 provides arcsecond-level accuracy for dates within a few thousand
/// years of J2000.0. For higher precision or dates far from J2000.0, consider
/// using JPL Development Ephemerides.
///
/// ## Usage
///
/// ```rust
/// use siderust::coordinates::transform::ephemeris::{BodyId, BodyEphemeris, Vsop87Ephemeris};
/// use siderust::astro::JulianDate;
///
/// let eph = Vsop87Ephemeris;
/// let earth_pos = eph.position_barycentric(BodyId::Earth, JulianDate::J2000);
/// println!("Earth position: {:?} AU", earth_pos);
/// ```
///
/// ## Zero-Cost Abstraction
///
/// `Vsop87Ephemeris` is a zero-sized type (unit struct). Using it has no runtime
/// overhead compared to calling VSOP87 functions directly.
#[derive(Debug, Clone, Copy, Default)]
pub struct Vsop87Ephemeris;

impl BodyEphemeris for Vsop87Ephemeris {
    fn position_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        match body {
            BodyId::Sun => {
                let target = Sun::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Earth => {
                let target = Earth::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Mercury => {
                let target = Mercury::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Venus => {
                let target = Venus::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Mars => {
                let target = Mars::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Jupiter => {
                let target = Jupiter::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Saturn => {
                let target = Saturn::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Uranus => {
                let target = Uranus::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Neptune => {
                let target = Neptune::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Moon => {
                // Moon is not in VSOP87; use ELP2000 via Earth-Moon barycenter
                // For now, approximate as Earth position (Moon orbit is small)
                // TODO: Implement proper lunar ephemeris
                let target = Earth::vsop87e(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
        }
    }

    fn position_heliocentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        match body {
            BodyId::Sun => [0.0, 0.0, 0.0],
            BodyId::Earth => {
                let target = Earth::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Mercury => {
                let target = Mercury::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Venus => {
                let target = Venus::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Mars => {
                let target = Mars::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Jupiter => {
                let target = Jupiter::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Saturn => {
                let target = Saturn::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Uranus => {
                let target = Uranus::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Neptune => {
                let target = Neptune::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
            BodyId::Moon => {
                // Approximate Moon as at Earth's position for now
                let target = Earth::vsop87a(jd);
                let pos = target.get_position();
                [pos.x().value(), pos.y().value(), pos.z().value()]
            }
        }
    }
}

impl VelocityEphemeris for Vsop87Ephemeris {
    fn velocity_barycentric(&self, body: BodyId, jd: JulianDate) -> [f64; 3] {
        match body {
            BodyId::Earth => {
                let vel = Earth::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Mercury => {
                let vel = Mercury::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Venus => {
                let vel = Venus::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Mars => {
                let vel = Mars::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Jupiter => {
                let vel = Jupiter::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Saturn => {
                let vel = Saturn::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Uranus => {
                let vel = Uranus::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            BodyId::Neptune => {
                let vel = Neptune::vsop87e_vel(jd);
                [vel.x().value(), vel.y().value(), vel.z().value()]
            }
            // Sun has negligible velocity in barycentric frame
            BodyId::Sun | BodyId::Moon => [0.0, 0.0, 0.0],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vsop87_earth_position() {
        let eph = Vsop87Ephemeris;
        let pos = eph.position_barycentric(BodyId::Earth, JulianDate::J2000);

        // Earth should be roughly 1 AU from barycenter
        let dist = (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt();
        assert!((dist - 1.0).abs() < 0.02, "Earth ~1 AU from barycenter");
    }

    #[test]
    fn test_vsop87_sun_heliocentric_is_zero() {
        let eph = Vsop87Ephemeris;
        let pos = eph.position_heliocentric(BodyId::Sun, JulianDate::J2000);

        assert_eq!(pos, [0.0, 0.0, 0.0], "Sun at heliocentric origin");
    }

    #[test]
    fn test_vsop87_earth_velocity() {
        let eph = Vsop87Ephemeris;
        let vel = eph.velocity_barycentric(BodyId::Earth, JulianDate::J2000);

        // Earth orbital velocity is roughly 0.017 AU/day
        let speed = (vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2)).sqrt();
        assert!(
            (speed - 0.017).abs() < 0.002,
            "Earth orbital speed ~0.017 AU/day, got {}",
            speed
        );
    }

    #[test]
    fn test_consistency_with_direct_calls() {
        let eph = Vsop87Ephemeris;
        let jd = JulianDate::J2000;

        // Compare trait method with direct VSOP87 call
        let via_trait = eph.position_barycentric(BodyId::Earth, jd);
        let target = Earth::vsop87e(jd);
        let direct = target.get_position();

        assert!((via_trait[0] - direct.x().value()).abs() < 1e-15);
        assert!((via_trait[1] - direct.y().value()).abs() < 1e-15);
        assert!((via_trait[2] - direct.z().value()).abs() < 1e-15);
    }
}
