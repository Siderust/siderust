// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Aberration (observer-velocity)
//!
//! This module applies the special-relativistic aberration of light due to an
//! observer's velocity.
//!
//! ```text
//! max effect   ≃ 20.5″
//! annual (Earth orbital) ≃ 20.5″
//! diurnal (Earth rotation) ≃ 0.3″
//! ```
//!
//! ## Velocity model
//! * The convenience functions in this module use **VSOP87E barycentric Earth
//!   velocity** (SSB-referenced), rotated from the dynamical *ecliptic J2000*
//!   frame into [`frames::EquatorialMeanJ2000`].
//!
//! ## References
//! * IERS Conventions 2020, §7.2 (aberration)
//! * SOFA/ERFA: stellar aberration is a full SR (Lorentz) transform (see e.g.
//!   `iauAb` / `eraAb` behavior)
//!
//! ## Implementation notes
//! * Uses the full special-relativistic aberration formula (Lorentz transform).
//! * Uses exact SI definitions for `c`, day, and AU to compute `c` in AU/day.

use crate::bodies::solar_system::Earth;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{
    cartesian::{direction, position, Velocity},
    frames,
};
use crate::time::JulianDate;
use qtty::*;

type AuPerDay = qtty::Per<AstronomicalUnit, Day>;
type AusPerDay = qtty::velocity::Velocity<AstronomicalUnit, Day>;

/// Speed of light in AU/day from exact SI definitions:
/// `c = 299_792_458 m/s`, `day = 86_400 s`, `AU = 149_597_870_700 m`.
pub const AU_PER_DAY_C_F64: f64 = 173.144_632_674_240_33_f64;

/// Same as [`AU_PER_DAY_C_F64`], as a quantity.
pub const AU_PER_DAY_C: AusPerDay = AusPerDay::new(AU_PER_DAY_C_F64);

#[inline]
fn aberrate_unit_vector_lorentz(
    u: nalgebra::Vector3<f64>,
    beta: nalgebra::Vector3<f64>,
) -> nalgebra::Vector3<f64> {
    let beta2 = beta.dot(&beta);
    if beta2 == 0.0 {
        return u;
    }

    // gamma = 1 / sqrt(1 - |beta|^2)
    let gamma = 1.0 / (1.0 - beta2).sqrt();
    let beta_dot_u = beta.dot(&u);

    // u' = [ u/gamma + beta + (gamma/(gamma+1)) (beta·u) beta ] / (1 + beta·u)
    let factor = gamma / (gamma + 1.0);
    let numerator = (u / gamma) + beta * (1.0 + factor * beta_dot_u);
    let denom = 1.0 + beta_dot_u;
    numerator / denom
}

/// Apply aberration to a unit direction in [`frames::EquatorialMeanJ2000`] using an explicit observer velocity.
#[must_use]
pub fn apply_aberration_to_direction_with_velocity(
    mean: direction::EquatorialMeanJ2000,
    velocity: &Velocity<frames::EquatorialMeanJ2000, AuPerDay>,
) -> direction::EquatorialMeanJ2000 {
    let beta = nalgebra::Vector3::new(
        velocity.x().value() / AU_PER_DAY_C_F64,
        velocity.y().value() / AU_PER_DAY_C_F64,
        velocity.z().value() / AU_PER_DAY_C_F64,
    );
    let u = nalgebra::Vector3::new(mean.x(), mean.y(), mean.z());
    let up = aberrate_unit_vector_lorentz(u, beta);
    direction::EquatorialMeanJ2000::normalize(up.x, up.y, up.z)
}

/// Remove aberration from a unit direction in [`frames::EquatorialMeanJ2000`] using an explicit observer velocity.
#[must_use]
pub fn remove_aberration_from_direction_with_velocity(
    app: direction::EquatorialMeanJ2000,
    velocity: &Velocity<frames::EquatorialMeanJ2000, AuPerDay>,
) -> direction::EquatorialMeanJ2000 {
    // Inverse is the same Lorentz transform with negated velocity.
    let beta = nalgebra::Vector3::new(
        -velocity.x().value() / AU_PER_DAY_C_F64,
        -velocity.y().value() / AU_PER_DAY_C_F64,
        -velocity.z().value() / AU_PER_DAY_C_F64,
    );
    let u = nalgebra::Vector3::new(app.x(), app.y(), app.z());
    let up = aberrate_unit_vector_lorentz(u, beta);
    direction::EquatorialMeanJ2000::normalize(up.x, up.y, up.z)
}

/// Apply **annual aberration** to a unit direction vector (mean J2000).
///
/// * `mean` – Geocentric unit vector in the mean equator & equinox of J2000.
/// * `jd`   – Epoch for Earth state evaluation (VSOP87 expects TDB; TT is a close approximation).
///
/// Returns a new [`crate::coordinates::cartesian::direction::EquatorialMeanJ2000`] including annual aberration.
#[must_use]
pub fn apply_aberration_to_direction(
    mean: direction::EquatorialMeanJ2000,
    jd: JulianDate,
) -> direction::EquatorialMeanJ2000 {
    // Use SSB-referenced (barycentric) Earth velocity for annual aberration.
    let velocity_ecl = Earth::vsop87e_vel(jd);
    let velocity: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = velocity_ecl.to_frame();
    apply_aberration_to_direction_with_velocity(mean, &velocity)
}

/// Remove **annual aberration** from an apparent direction.
/// Inverse operation of [`apply_aberration_to_direction`].
#[must_use]
pub fn remove_aberration_from_direction(
    app: direction::EquatorialMeanJ2000,
    jd: JulianDate,
) -> direction::EquatorialMeanJ2000 {
    let velocity_ecl = Earth::vsop87e_vel(jd);
    let velocity: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = velocity_ecl.to_frame();
    remove_aberration_from_direction_with_velocity(app, &velocity)
}

/// Apply **annual aberration** to a position vector, preserving its
/// geocentric distance.
#[must_use]
pub fn apply_aberration<U: LengthUnit>(
    mean: position::EquatorialMeanJ2000<U>,
    jd: JulianDate,
) -> position::EquatorialMeanJ2000<U> {
    if mean.distance() == 0.0 {
        // Don't look at your feet!
        return mean;
    }

    // Safe to unwrap: we just checked distance is non-zero
    let dir = mean
        .direction()
        .expect("non-zero position should have a direction");
    apply_aberration_to_direction(dir, jd).position(mean.distance())
}

/// Remove **annual aberration** from a position vector, preserving its
/// geocentric distance.
#[must_use]
pub fn remove_aberration<U: LengthUnit>(
    app: position::EquatorialMeanJ2000<U>,
    jd: JulianDate,
) -> position::EquatorialMeanJ2000<U> {
    if app.distance() == 0.0 {
        // Don't look at your feet!
        return app;
    }

    // Safe to unwrap: we just checked distance is non-zero
    let dir = app
        .direction()
        .expect("non-zero position should have a direction");
    remove_aberration_from_direction(dir, jd).position(app.distance())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::spherical::{self, position};
    use nalgebra::Vector3;

    fn apply_aberration_sph<U: LengthUnit>(
        mean: &position::EquatorialMeanJ2000<U>,
        jd: JulianDate,
    ) -> position::EquatorialMeanJ2000<U> {
        spherical::Position::from_cartesian(&apply_aberration(mean.to_cartesian(), jd))
    }

    #[test]
    fn test_aberration_preserva_distance_and_epoch() {
        let jd = JulianDate::new(2451545.0); // J2000.0
        let mean =
            position::EquatorialMeanJ2000::<Au>::new(Degrees::new(10.0), Degrees::new(20.0), 1.23);
        let out = apply_aberration_sph(&mean, jd);

        assert_eq!(out.distance(), mean.distance());
    }

    #[test]
    fn test_aberration_introduces_shift() {
        let jd = JulianDate::new(2451545.0); // J2000.0
        let mean = position::EquatorialMeanJ2000::<Au>::new(
            Degrees::new(0.0), // RA = 0°
            Degrees::new(0.0), // Dec = 0°
            1.0,
        );
        let out = apply_aberration_sph(&mean, jd);

        let delta_ra = out.ra().abs_separation(mean.ra());
        let delta_dec = out.dec().abs_separation(mean.dec());
        assert!(
            delta_ra > 0.0 || delta_dec > 0.0,
            "Expected a change in RA or Dec"
        );
        assert!(delta_ra < 0.01 && delta_dec < 0.01, "Shift is too large")
    }

    #[test]
    fn test_aberration_at_north_pole() {
        let jd = JulianDate::new(2451545.0);
        let mean = position::EquatorialMeanJ2000::<Au>::new(
            Degrees::new(123.4), // dummy RA
            Degrees::new(90.0),  // Dec = +90°
            1.0,
        );
        let out = apply_aberration_sph(&mean, jd);

        assert!(
            out.dec() < 90.0,
            "Declination should decrease slightly at pole"
        );
        assert!(!out.ra().is_nan(), "RA must not be NaN at the pole");
    }

    #[test]
    fn test_speed_of_light() {
        // Exact from SI definitions (to ~1e-15 relative precision in f64):
        // 299792458 * 86400 / 149597870700 = 173.14463267424033...
        let expected = AusPerDay::new(173.144_632_674_240_33);
        assert!((AU_PER_DAY_C - expected).abs() < AusPerDay::new(1e-12));
    }

    #[test]
    fn aberration_roundtrip_is_machine_precision() {
        let jd = JulianDate::new(2458850.0); // 2020-ish
        let velocity_ecl = Earth::vsop87e_vel(jd);
        let velocity: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = velocity_ecl.to_frame();

        let directions = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.3, -0.7, 0.64),
            Vector3::new(-0.8, 0.2, -0.56),
        ];

        for u in directions {
            let mean = direction::EquatorialMeanJ2000::from_vec3(u.normalize());
            let app = apply_aberration_to_direction_with_velocity(mean, &velocity);
            let rec = remove_aberration_from_direction_with_velocity(app, &velocity);

            let dot = mean.x() * rec.x() + mean.y() * rec.y() + mean.z() * rec.z();
            let ang = dot.clamp(-1.0, 1.0).acos();
            assert!(ang < 5e-15, "roundtrip angle too large: {}", ang);
        }
    }
}
