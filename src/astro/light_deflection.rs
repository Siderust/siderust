// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Gravitational Light Deflection
//!
//! When light from a distant source passes near a massive body (primarily the
//! Sun), its path is curved by the gravitational field. This shifts the
//! **apparent** direction of the source away from the deflecting body.
//!
//! ## Solar deflection
//!
//! The maximum solar deflection at the limb is **1.75″** (the classic Einstein
//! value). For a source at angular distance **θ** from the Sun:
//!
//! ```text
//! Δθ ≈ (1.75″ / sin θ) × (1 + cos θ) / 2
//! ```
//!
//! or equivalently in the vector formulation (IERS Conventions 2010, §7.1.1):
//!
//! ```text
//! δs = (2 G M / c² |q|) × [ (s · q̂) q̂ − (q̂ · s) s ] / (1 + q̂ · s)
//! ```
//!
//! where:
//! - **s**: unit direction vector toward the source (BCRS)
//! - **q**: vector from the deflecting body to the observer (BCRS)
//! - **q̂** = q / |q|
//!
//! ## Effect magnitudes
//!
//! | Body   | Max deflection |
//! |--------|---------------|
//! | Sun    | 1.75″ (limb)  |
//! | Jupiter| 0.017″        |
//! | Saturn | 0.006″        |
//! | Moon   | 0.026 mas     |
//!
//! For most applications, only the Sun's contribution matters above 1 mas.
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.6
//! * IERS Conventions (2010), §7.1.1
//! * SOFA routine `iauLdsun`, `iauLd`
//! * Klioner, S. A. (2003), AJ 125, 1580

use crate::coordinates::cartesian::direction;
use qtty::{Arcseconds, AstronomicalUnits, Radians};

/// Schwarzschild radius of the Sun in AU: 2 G M☉ / c².
///
/// Value: 1.97412574336e-8 AU (≈ 2953.25 m).
///
/// ## References
/// * IAU 2015 Resolution B3 (nominal solar mass parameter)
/// * IERS Conventions (2010), Table 1.1
const SOLAR_SCHWARZSCHILD_AU: f64 = 1.974_125_743_36e-8;

/// Apply solar gravitational light deflection to a source direction.
///
/// Given the unit direction vector of a source and the Sun-observer geometry,
/// this computes the apparent direction after accounting for the Sun's
/// gravitational field bending light rays.
///
/// ## Parameters
///
/// * `star`: Unit direction toward the source in [`EquatorialMeanJ2000`](direction::EquatorialMeanJ2000) (BCRS).
/// * `earth_sun_au`: Vector from the Sun to the Earth/observer in AU `[x, y, z]` (BCRS).
///
/// ## Returns
///
/// The deflected unit direction in [`EquatorialMeanJ2000`](direction::EquatorialMeanJ2000).
///
/// ## Notes
///
/// - Sources within ~5° of the Sun may need higher-order corrections not
///   included here.
/// - The deflection is computed in the BCRS; for GCRS applications the
///   difference is negligible (< 1 μas).
///
/// ## References
/// * SOFA routine `iauLdsun`
/// * IERS Conventions (2010), eq. 7.4
pub fn solar_deflection(
    star: direction::EquatorialMeanJ2000,
    earth_sun_au: [f64; 3],
) -> direction::EquatorialMeanJ2000 {
    // ── Geometry ──
    let q_mag =
        (earth_sun_au[0].powi(2) + earth_sun_au[1].powi(2) + earth_sun_au[2].powi(2)).sqrt();

    if q_mag < 1e-10 {
        return star; // degenerate: observer at the Sun
    }

    // Unit vector from Sun toward observer
    let q_hat = [
        earth_sun_au[0] / q_mag,
        earth_sun_au[1] / q_mag,
        earth_sun_au[2] / q_mag,
    ];

    // Dot products
    let sx = star.x();
    let sy = star.y();
    let sz = star.z();
    let s_dot_q = sx * q_hat[0] + sy * q_hat[1] + sz * q_hat[2];

    // Avoid degenerate case: source is exactly behind the Sun
    let denom = 1.0 + s_dot_q;
    if denom.abs() < 1e-15 {
        return star; // source directly behind sun — deflection undefined
    }

    // Deflection factor: 2 G M / (c² |q|)
    let factor = SOLAR_SCHWARZSCHILD_AU / q_mag;

    // Vector deflection formula (SOFA iauLd):
    // p1 = p + f * (qhat - sn*p) where f = 2mu/(c^2 * eq) / (1 + sn)
    // and sn = p · qhat
    let f = factor / denom;

    let dx = sx + f * (q_hat[0] - s_dot_q * sx);
    let dy = sy + f * (q_hat[1] - s_dot_q * sy);
    let dz = sz + f * (q_hat[2] - s_dot_q * sz);

    // Re-normalize to unit direction
    direction::EquatorialMeanJ2000::normalize(dx, dy, dz)
}

/// Compute the magnitude of solar deflection for a source at a given
/// angular distance from the Sun.
///
/// ## Parameters
///
/// * `sun_angle`: Angular distance of the source from the Sun.
/// * `sun_distance`: Observer-Sun distance (typically ~1 AU).
///
/// ## Returns
///
/// The deflection angle in [`Arcsecond`].
///
/// ## Notes
///
/// Uses the classical formula: Δθ = 1.75″ × (1/tan(θ/2)) × (1/R)
/// where R is the Sun distance in AU.
#[inline]
pub fn solar_deflection_magnitude(
    sun_angle: Radians,
    sun_distance: AstronomicalUnits,
) -> Arcseconds {
    // Maximum deflection at limb at 1 AU: 1.7512 arcsec
    // (From 2 G M☉ / (c² R☉_angle))
    const LIMB_DEFLECTION_1AU: f64 = 1.7512;

    let angle = sun_angle.value();
    let dist = sun_distance.value();

    if angle.abs() < 1e-10 {
        return Arcseconds::new(0.0); // directly at the Sun — meaningless
    }

    // Δθ = (1.75″/R) × cos(θ/2) / sin(θ/2)  [cotangent form]
    // This avoids singularity at θ = 0 better than 1/sinθ form
    let half = angle / 2.0;
    Arcseconds::new(LIMB_DEFLECTION_1AU / dist * half.cos() / half.sin().max(1e-10))
}

/// Remove solar deflection from an apparent direction to get the geometric
/// (catalog) direction.
///
/// This is the inverse of [`solar_deflection`]: given the apparent (observed)
/// direction, compute the un-deflected direction.
///
/// For small deflections (which is always the case for the Sun at > a few
/// degrees), one iteration suffices. For sources very near the limb,
/// a second iteration can be applied.
///
/// ## Parameters
///
/// * `apparent`: Apparent unit direction in [`EquatorialMeanJ2000`](direction::EquatorialMeanJ2000) (BCRS).
/// * `earth_sun_au`: Vector from the Sun to the Earth/observer in AU `[x, y, z]` (BCRS).
///
/// ## Returns
///
/// The un-deflected unit direction in [`EquatorialMeanJ2000`](direction::EquatorialMeanJ2000).
pub fn solar_deflection_inverse(
    apparent: direction::EquatorialMeanJ2000,
    earth_sun_au: [f64; 3],
) -> direction::EquatorialMeanJ2000 {
    // Single iteration: apply negative deflection
    // For small angles: geometric ≈ apparent − δ(apparent)
    let deflected = solar_deflection(apparent, earth_sun_au);

    let gx = 2.0 * apparent.x() - deflected.x();
    let gy = 2.0 * apparent.y() - deflected.y();
    let gz = 2.0 * apparent.z() - deflected.z();

    // Re-normalize
    direction::EquatorialMeanJ2000::normalize(gx, gy, gz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Radians;

    #[test]
    fn deflection_at_limb() {
        // Source 90° from Sun, observer at 1 AU
        // Deflection should be ≈ 1.75″ / tan(45°) ≈ 1.75″
        let angle = Radians::new(std::f64::consts::FRAC_PI_2);
        let defl = solar_deflection_magnitude(angle, AstronomicalUnits::new(1.0));
        assert!(
            (defl.value() - 1.75).abs() < 0.1,
            "deflection at 90° = {}″, expected ≈ 1.75″",
            defl.value()
        );
    }

    #[test]
    fn deflection_far_from_sun() {
        // Source 90° away from the Sun
        let angle = Radians::new(std::f64::consts::FRAC_PI_2);
        let defl = solar_deflection_magnitude(angle, AstronomicalUnits::new(1.0));
        // At 90°, cot(45°) = 1.0, so deflection ≈ 1.75″
        assert!(defl.value() > 0.0 && defl.value() < 10.0);
    }

    #[test]
    fn deflection_decreases_with_distance() {
        let angle = Radians::new(0.5); // ~29°
        let d1 = solar_deflection_magnitude(angle, AstronomicalUnits::new(1.0));
        let d2 = solar_deflection_magnitude(angle, AstronomicalUnits::new(2.0));
        assert!(
            d1.value() > d2.value(),
            "deflection should decrease with Sun distance"
        );
        assert!(
            (d1.value() / d2.value() - 2.0).abs() < 0.01,
            "should scale as 1/R"
        );
    }

    #[test]
    fn vector_deflection_direction() {
        // Source at (1, 0, 0), Sun along (0, 0, -1) direction from observer
        // observer at 1 AU from Sun
        let star = direction::EquatorialMeanJ2000::new(1.0, 0.0, 0.0);
        let earth_sun_au = [0.0, 0.0, -1.0]; // Sun 1 AU in -z direction

        let deflected = solar_deflection(star, earth_sun_au);

        // Deflection should be tiny: ~1.75″ ≈ 8.5e-6 rad
        // The source is 90° from the Sun, so deflection ≈ 1.75″
        let angular_diff = (deflected.x() - star.x()).powi(2)
            + (deflected.y() - star.y()).powi(2)
            + (deflected.z() - star.z()).powi(2);
        let angular_diff = angular_diff.sqrt();
        // Should be on the order of microradians
        assert!(
            angular_diff > 1e-9 && angular_diff < 1e-4,
            "deflection magnitude = {} rad, expected ~8.5e-6",
            angular_diff
        );
    }

    #[test]
    fn no_deflection_when_at_sun() {
        let star = direction::EquatorialMeanJ2000::new(1.0, 0.0, 0.0);
        let earth_sun_au = [0.0, 0.0, 0.0]; // degenerate
        let result = solar_deflection(star, earth_sun_au);
        assert_eq!(result.x(), star.x());
        assert_eq!(result.y(), star.y());
        assert_eq!(result.z(), star.z());
    }

    #[test]
    fn deflection_roundtrip() {
        // Apply deflection, then inverse — should recover original direction
        let star = direction::EquatorialMeanJ2000::new(0.6, 0.7, 0.3742);
        let earth_sun_au = [0.5, -0.3, 0.8]; // ~1 AU

        let deflected = solar_deflection(star, earth_sun_au);
        let recovered = solar_deflection_inverse(deflected, earth_sun_au);

        assert!(
            (recovered.x() - star.x()).abs() < 1e-10,
            "roundtrip mismatch x: {} vs {}",
            recovered.x(),
            star.x()
        );
        assert!(
            (recovered.y() - star.y()).abs() < 1e-10,
            "roundtrip mismatch y: {} vs {}",
            recovered.y(),
            star.y()
        );
        assert!(
            (recovered.z() - star.z()).abs() < 1e-10,
            "roundtrip mismatch z: {} vs {}",
            recovered.z(),
            star.z()
        );
    }
}
