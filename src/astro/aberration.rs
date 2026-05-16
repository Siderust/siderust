// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Aberration (observer-velocity)
//!
//! Applies the special-relativistic aberration of light to celestial
//! directions, given an observer's velocity in a barycentric frame.
//!
//! ## Scientific scope
//!
//! Stellar aberration is the apparent shift of a source's direction caused by
//! the finite speed of light combined with the observer's motion relative to
//! the rest frame of the source. The annual (Earth orbital) component reaches
//! ≈ 20.5″, while the diurnal (Earth rotation) component is ≈ 0.3″. Correctly
//! modelling aberration is essential for sub-arcsecond astrometry and for any
//! observer-relative direction that is to be compared with catalogue (BCRS)
//! positions.
//!
//! ## Technical scope
//!
//! The convenience functions in this module use the **VSOP87E barycentric
//! Earth velocity** (SSB-referenced), rotated from the dynamical *ecliptic
//! J2000* frame into [`frames::EquatorialMeanJ2000`]. The aberration is
//! applied via the full Lorentz transform of unit direction vectors rather
//! than the classical first-order approximation, ensuring < 1 μas residuals.
//! The speed of light in AU/day is computed from exact SI definitions of
//! `c`, the day, and the AU.
//!
//! ## References
//!
//! * IERS Conventions 2020, §7.2 (aberration)
//! * SOFA: stellar aberration as a full SR (Lorentz) transform
//!   (cf. `iauAb` / `eraAb`)

use crate::calculus::ephemeris::Ephemeris;
use crate::coordinates::transform::context::DefaultEphemeris;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{
    cartesian::{direction, position, Velocity},
    frames,
};
use crate::qtty::*;
use crate::time::JulianDate;

type AuPerDay = crate::qtty::Per<AstronomicalUnit, Day>;

/// Full SR Lorentz aberration transform for a unit direction vector.
///
/// `beta` is the observer's velocity divided by c (dimensionless, non-unit-length).
/// The Lorentz formula is:
/// `u' = [ u/γ + β (1 + γ/(γ+1) β·u) ] / (1 + β·u)`
#[inline]
fn aberrate_direction_lorentz<F: frames::ReferenceFrame>(
    u: direction::Direction<F>,
    beta: Velocity<F, Ratio>,
) -> direction::Direction<F> {
    let beta_raw = affn::cartesian::XYZ::from_array(*beta.as_array()).to_raw();
    let beta2 = beta_raw.magnitude_squared();
    if beta2 == 0.0 {
        return u;
    }

    let gamma = 1.0 / (1.0 - beta2).sqrt();
    let u_vec = Velocity::<F, Ratio>::new(u.x(), u.y(), u.z());
    let beta_dot_u = beta_raw.dot(&affn::cartesian::XYZ::from_array(u.as_array()));

    let result = (u_vec.scale(1.0 / gamma)
        + beta.scale(1.0 + (gamma / (gamma + 1.0)) * beta_dot_u))
        / Quantity::<Ratio>::new(1.0 + beta_dot_u);

    direction::Direction::from_array(
        affn::cartesian::XYZ::from_array(*result.as_array())
            .to_raw()
            .into_array(),
    )
}

/// Apply aberration to a unit direction in [`frames::EquatorialMeanJ2000`] using an explicit observer velocity.
#[must_use]
pub fn apply_aberration_to_direction_with_velocity(
    mean: direction::EquatorialMeanJ2000,
    velocity: &Velocity<frames::EquatorialMeanJ2000, AuPerDay>,
) -> direction::EquatorialMeanJ2000 {
    aberrate_direction_lorentz(mean, *velocity / crate::qtty::velocity::AU_PER_DAY_C)
}

/// Remove aberration from a unit direction in [`frames::EquatorialMeanJ2000`] using an explicit observer velocity.
#[must_use]
pub fn remove_aberration_from_direction_with_velocity(
    app: direction::EquatorialMeanJ2000,
    velocity: &Velocity<frames::EquatorialMeanJ2000, AuPerDay>,
) -> direction::EquatorialMeanJ2000 {
    // Inverse is the same Lorentz transform with negated velocity.
    aberrate_direction_lorentz(app, -(*velocity / crate::qtty::velocity::AU_PER_DAY_C))
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
    let velocity_ecl = DefaultEphemeris::earth_barycentric_velocity(jd);
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
    let velocity_ecl = DefaultEphemeris::earth_barycentric_velocity(jd);
    let velocity: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = velocity_ecl.to_frame();
    remove_aberration_from_direction_with_velocity(app, &velocity)
}

/// Apply **annual aberration** to a position vector, preserving its
/// geocentric distance.
///
/// # Panics
///
/// Will not panic on user input: the zero-distance case (`|p| == 0`) is
/// short-circuited early to return the input unchanged. The internal
/// `direction()` extraction is guarded by that distance check, so the
/// `expect` it contains is unreachable for any well-formed `Position`.
#[must_use]
pub fn apply_aberration<U: LengthUnit>(
    mean: position::EquatorialMeanJ2000<U>,
    jd: JulianDate,
) -> position::EquatorialMeanJ2000<U> {
    if mean.distance() == Quantity::<U>::new(0.0) {
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
///
/// # Panics
///
/// As for [`apply_aberration`]: the zero-distance case is short-circuited
/// and the internal direction extraction is provably safe under the
/// guarded distance check.
#[must_use]
pub fn remove_aberration<U: LengthUnit>(
    app: position::EquatorialMeanJ2000<U>,
    jd: JulianDate,
) -> position::EquatorialMeanJ2000<U> {
    if app.distance() == Quantity::<U>::new(0.0) {
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

    type AusPerDay = crate::qtty::velocity::Velocity<AstronomicalUnit, Day>;

    fn apply_aberration_sph<U: LengthUnit>(
        mean: &position::EquatorialMeanJ2000<U>,
        jd: JulianDate,
    ) -> position::EquatorialMeanJ2000<U> {
        spherical::Position::from_cartesian(&apply_aberration(mean.to_cartesian(), jd))
    }

    #[test]
    fn test_aberration_preserva_distance_and_epoch() {
        let jd = crate::time::jd(qtty::Day::new(2451545.0)); // J2000.0
        let mean =
            position::EquatorialMeanJ2000::<Au>::new(Degrees::new(10.0), Degrees::new(20.0), 1.23);
        let out = apply_aberration_sph(&mean, jd);

        // Distance should be preserved to within floating-point rounding
        assert!(
            (out.distance - mean.distance).abs() < AstronomicalUnits::new(1e-12),
            "Distance should be preserved: got {:?}, expected {:?}",
            out.distance,
            mean.distance
        );
    }

    #[test]
    fn test_aberration_introduces_shift() {
        let jd = crate::time::jd(qtty::Day::new(2451545.0)); // J2000.0
        let mean = position::EquatorialMeanJ2000::<Au>::new(
            Degrees::new(0.0), // RA = 0°
            Degrees::new(0.0), // Dec = 0°
            1.0,
        );
        let out = apply_aberration_sph(&mean, jd);

        let delta_ra = out.ra().abs_separation(mean.ra());
        let delta_dec = out.dec().abs_separation(mean.dec());
        assert!(
            delta_ra > Degrees::new(0.0) || delta_dec > Degrees::new(0.0),
            "Expected a change in RA or Dec"
        );
        assert!(
            delta_ra < Degrees::new(0.01) && delta_dec < Degrees::new(0.01),
            "Shift is too large"
        )
    }

    #[test]
    fn test_aberration_at_north_pole() {
        let jd = crate::time::jd(qtty::Day::new(2451545.0));
        let mean = position::EquatorialMeanJ2000::<Au>::new(
            Degrees::new(123.4), // dummy RA
            Degrees::new(90.0),  // Dec = +90°
            1.0,
        );
        let out = apply_aberration_sph(&mean, jd);

        assert!(
            out.dec() < Degrees::new(90.0),
            "Declination should decrease slightly at pole"
        );
        assert!(!out.ra().is_nan(), "RA must not be NaN at the pole");
    }

    #[test]
    fn test_speed_of_light() {
        // Exact from SI definitions (to ~1e-15 relative precision in f64):
        // 299792458 * 86400 / 149597870700 = 173.14463267424033...
        let expected = AusPerDay::new(173.144_632_674_240_33);
        assert!((crate::qtty::velocity::AU_PER_DAY_C - expected).abs() < AusPerDay::new(1e-12));
    }

    #[test]
    fn aberration_roundtrip_is_machine_precision() {
        use crate::bodies::solar_system::Earth;
        let jd = crate::time::jd(qtty::Day::new(2458850.0)); // 2020-ish
        let velocity_ecl = Earth::vsop87e_vel(jd);
        let velocity: Velocity<frames::EquatorialMeanJ2000, AuPerDay> = velocity_ecl.to_frame();

        let directions: [[f64; 3]; 5] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.3, -0.7, 0.64],
            [-0.8, 0.2, -0.56],
        ];

        for u in directions {
            let mean = direction::EquatorialMeanJ2000::from_array(u);
            let app = apply_aberration_to_direction_with_velocity(mean, &velocity);
            let rec = remove_aberration_from_direction_with_velocity(app, &velocity);

            let dot = mean.x() * rec.x() + mean.y() * rec.y() + mean.z() * rec.z();
            let ang = dot.clamp(-1.0, 1.0).acos();
            assert!(ang < 5e-15, "roundtrip angle too large: {}", ang);
        }
    }
}
