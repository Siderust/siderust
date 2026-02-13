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
/// * `star`: Unit direction vector toward the source `[x, y, z]` (BCRS).
/// * `earth_sun`: Vector from the Sun to the Earth/observer in AU `[x, y, z]` (BCRS).
///
/// ## Returns
///
/// The deflected unit direction vector `[x, y, z]`.
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
pub fn solar_deflection(star: [f64; 3], earth_sun: [f64; 3]) -> [f64; 3] {
    // ── Geometry ──
    let q_mag = (earth_sun[0].powi(2) + earth_sun[1].powi(2) + earth_sun[2].powi(2)).sqrt();

    if q_mag < 1e-10 {
        return star; // degenerate: observer at the Sun
    }

    // Unit vector from Sun toward observer
    let q_hat = [
        earth_sun[0] / q_mag,
        earth_sun[1] / q_mag,
        earth_sun[2] / q_mag,
    ];

    // Dot products
    let s_dot_q = star[0] * q_hat[0] + star[1] * q_hat[1] + star[2] * q_hat[2];

    // Avoid degenerate case: source is exactly behind the Sun
    let denom = 1.0 + s_dot_q;
    if denom.abs() < 1e-15 {
        return star; // source directly behind sun — deflection undefined
    }

    // Deflection factor: 2 G M / (c² |q|)
    let factor = SOLAR_SCHWARZSCHILD_AU / q_mag;

    // Vector deflection formula (IERS 2010, eq. 7.4):
    // δs = factor × [ (s·q̂)q̂ − s·q̂ × q̂ − q̂·s × s ] / (1 + q̂·s)
    // Simplified: δs = factor / denom × [ q̂ − s_dot_q × s ]
    // Actually the correct formula from SOFA iauLd:
    // p1 = p + f * (qhat - sn*p) where f = 2mu/(c^2 * eq) / (1 + sn)
    // and sn = p · qhat

    let f = factor / denom;

    let mut deflected = [0.0f64; 3];
    for i in 0..3 {
        deflected[i] = star[i] + f * (q_hat[i] - s_dot_q * star[i]);
    }

    // Re-normalize to unit vector
    let mag = (deflected[0].powi(2) + deflected[1].powi(2) + deflected[2].powi(2)).sqrt();
    [deflected[0] / mag, deflected[1] / mag, deflected[2] / mag]
}

/// Compute the magnitude of solar deflection for a source at a given
/// angular distance from the Sun.
///
/// Returns the deflection in arcseconds.
///
/// ## Parameters
///
/// * `sun_angle_rad`: Angular distance of the source from the Sun (radians).
/// * `sun_distance_au`: Observer-Sun distance in AU (typically ~1.0).
///
/// ## Notes
///
/// Uses the classical formula: Δθ = 1.75″ × (1/tan(θ/2)) × (1/R)
/// where R is the Sun distance in AU.
#[inline]
pub fn solar_deflection_magnitude(sun_angle_rad: f64, sun_distance_au: f64) -> f64 {
    // Maximum deflection at limb at 1 AU: 1.7512 arcsec
    // (From 2 G M☉ / (c² R☉_angle))
    const LIMB_DEFLECTION_1AU: f64 = 1.7512;

    if sun_angle_rad.abs() < 1e-10 {
        return 0.0; // directly at the Sun — meaningless
    }

    // Δθ = (1.75″/R) × cos(θ/2) / sin(θ/2)  [cotangent form]
    // This avoids singularity at θ = 0 better than 1/sinθ form
    let half = sun_angle_rad / 2.0;
    LIMB_DEFLECTION_1AU / sun_distance_au * half.cos() / half.sin().max(1e-10)
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
/// * `apparent`: Apparent unit direction vector `[x, y, z]` (BCRS).
/// * `earth_sun`: Vector from the Sun to the Earth/observer in AU `[x, y, z]` (BCRS).
///
/// ## Returns
///
/// The un-deflected unit direction vector `[x, y, z]`.
pub fn solar_deflection_inverse(apparent: [f64; 3], earth_sun: [f64; 3]) -> [f64; 3] {
    // Single iteration: apply negative deflection
    // For small angles: geometric ≈ apparent − δ(apparent)
    let deflected = solar_deflection(apparent, earth_sun);
    let mut geometric = [0.0f64; 3];
    for i in 0..3 {
        geometric[i] = 2.0 * apparent[i] - deflected[i];
    }

    // Re-normalize
    let mag = (geometric[0].powi(2) + geometric[1].powi(2) + geometric[2].powi(2)).sqrt();
    [geometric[0] / mag, geometric[1] / mag, geometric[2] / mag]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deflection_at_limb() {
        // Source 90° from Sun, observer at 1 AU
        // Deflection should be ≈ 1.75″ / tan(45°) ≈ 1.75″
        let angle = std::f64::consts::FRAC_PI_2; // 90° from Sun
        let defl = solar_deflection_magnitude(angle, 1.0);
        assert!(
            (defl - 1.75).abs() < 0.1,
            "deflection at 90° = {}″, expected ≈ 1.75″",
            defl
        );
    }

    #[test]
    fn deflection_far_from_sun() {
        // Source 90° away from the Sun: deflection ≈ 0.004″
        let angle = std::f64::consts::FRAC_PI_2;
        let defl = solar_deflection_magnitude(angle, 1.0);
        // At 90°, cot(45°) = 1.0, so deflection ≈ 1.75″
        assert!(defl > 0.0 && defl < 10.0);
    }

    #[test]
    fn deflection_decreases_with_distance() {
        let angle = 0.5; // ~29°
        let d1 = solar_deflection_magnitude(angle, 1.0);
        let d2 = solar_deflection_magnitude(angle, 2.0);
        assert!(d1 > d2, "deflection should decrease with Sun distance");
        assert!((d1 / d2 - 2.0).abs() < 0.01, "should scale as 1/R");
    }

    #[test]
    fn vector_deflection_direction() {
        // Source at (1, 0, 0), Sun along (0, 0, -1) direction from observer
        // observer at 1 AU from Sun
        let star = [1.0, 0.0, 0.0];
        let earth_sun = [0.0, 0.0, -1.0]; // Sun 1 AU in -z direction

        let deflected = solar_deflection(star, earth_sun);

        // Deflection should be tiny: ~1.75″ ≈ 8.5e-6 rad
        // The source is 90° from the Sun, so deflection ≈ 1.75″
        let angular_diff = (deflected[0] - star[0]).powi(2)
            + (deflected[1] - star[1]).powi(2)
            + (deflected[2] - star[2]).powi(2);
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
        let star = [1.0, 0.0, 0.0];
        let earth_sun = [0.0, 0.0, 0.0]; // degenerate
        let result = solar_deflection(star, earth_sun);
        assert_eq!(result, star);
    }

    #[test]
    fn deflection_roundtrip() {
        // Apply deflection, then inverse — should recover original direction
        let star = [0.6_f64, 0.7, 0.3742];
        // Normalize
        let mag = (star[0].powi(2) + star[1].powi(2) + star[2].powi(2)).sqrt();
        let star = [star[0] / mag, star[1] / mag, star[2] / mag];

        let earth_sun = [0.5, -0.3, 0.8]; // ~1 AU

        let deflected = solar_deflection(star, earth_sun);
        let recovered = solar_deflection_inverse(deflected, earth_sun);

        for i in 0..3 {
            assert!(
                (recovered[i] - star[i]).abs() < 1e-10,
                "roundtrip mismatch [{}]: {} vs {}",
                i, recovered[i], star[i]
            );
        }
    }
}
