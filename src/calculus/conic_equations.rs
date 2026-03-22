// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Heliocentric propagation helpers for conic and mean-motion orbit models.
//!
//! All functions in this module assume a **heliocentric** context: the
//! gravitational parameter is the Gaussian gravitational constant
//! (`k = 0.01720209895 AU^{3/2} d^{-1}`), which implicitly encodes
//! `GM_☉` in the AU-day system. The returned positions are heliocentric
//! ecliptic J2000.
//!
//! To propagate orbits around other central bodies, use the generic
//! gravitational parameter approach in a future extension.

use crate::astro::conic::{
    map_validation_error, ConicError, ConicKind, ConicOrbit, MeanMotionOrbit,
};
use crate::astro::orbit::OrientationTrig;
use crate::calculus::kepler_equations::solve_keplers_equation;
use crate::coordinates::cartesian::position::EclipticMeanJ2000;
use crate::time::JulianDate;
use qtty::*;

/// Gaussian gravitational constant `k` in AU^{3/2} d^{-1}.
///
/// Encodes `√(GM_☉)` in the heliocentric AU-day system. Using this
/// constant ties all propagation in this module to a Sun-centered
/// context.
const GAUSSIAN_GRAVITATIONAL_CONSTANT: f64 = 0.01720209895;
const HYPERBOLIC_TOLERANCE: f64 = 1e-14;
const MAX_HYPERBOLIC_ITERS: usize = 100;

fn solve_hyperbolic_anomaly(mean_anomaly_radians: f64, eccentricity: f64) -> f64 {
    let mut hyperbolic_anomaly = (mean_anomaly_radians / eccentricity).asinh();
    for _ in 0..MAX_HYPERBOLIC_ITERS {
        let sinh_h = hyperbolic_anomaly.sinh();
        let cosh_h = hyperbolic_anomaly.cosh();
        let delta = (eccentricity * sinh_h - hyperbolic_anomaly - mean_anomaly_radians)
            / (eccentricity * cosh_h - 1.0);
        hyperbolic_anomaly -= delta;
        if delta.abs() < HYPERBOLIC_TOLERANCE {
            break;
        }
    }
    hyperbolic_anomaly
}

/// Like [`rotate_to_ecliptic`] but uses precomputed sin/cos values from
/// [`OrientationTrig`].
#[inline]
pub(crate) fn rotate_to_ecliptic_precomputed(
    radius_au: f64,
    trig: &OrientationTrig,
    true_anomaly_radians: f64,
) -> EclipticMeanJ2000<AstronomicalUnit> {
    // argument of latitude u = ω + ν
    let (sin_nu, cos_nu) = true_anomaly_radians.sin_cos();
    let sin_u = trig.sin_omega() * cos_nu + trig.cos_omega() * sin_nu;
    let cos_u = trig.cos_omega() * cos_nu - trig.sin_omega() * sin_nu;
    let x = radius_au * (trig.cos_node() * cos_u - trig.sin_node() * sin_u * trig.cos_i());
    let y = radius_au * (trig.sin_node() * cos_u + trig.cos_node() * sin_u * trig.cos_i());
    let z = radius_au * sin_u * trig.sin_i();
    EclipticMeanJ2000::new(
        AstronomicalUnits::new(x),
        AstronomicalUnits::new(y),
        AstronomicalUnits::new(z),
    )
}

#[inline]
pub(crate) fn rotate_to_ecliptic(
    radius_au: f64,
    inclination: Degrees,
    argument_of_periapsis: Degrees,
    longitude_of_ascending_node: Degrees,
    true_anomaly_radians: f64,
) -> EclipticMeanJ2000<AstronomicalUnit> {
    let (sin_i, cos_i) = inclination.to::<Radian>().value().sin_cos();
    let (sin_u, cos_u) =
        (argument_of_periapsis.to::<Radian>().value() + true_anomaly_radians).sin_cos();
    let (sin_node, cos_node) = longitude_of_ascending_node.to::<Radian>().value().sin_cos();
    let x = radius_au * (cos_node * cos_u - sin_node * sin_u * cos_i);
    let y = radius_au * (sin_node * cos_u + cos_node * sin_u * cos_i);
    let z = radius_au * sin_u * sin_i;
    EclipticMeanJ2000::new(
        AstronomicalUnits::new(x),
        AstronomicalUnits::new(y),
        AstronomicalUnits::new(z),
    )
}

/// Given a solved eccentric anomaly and eccentricity, returns the true anomaly
/// and the heliocentric radius for an elliptic orbit.
#[inline]
pub(crate) fn elliptic_true_anomaly_and_radius(
    eccentric_anomaly: Radians,
    eccentricity: f64,
    semi_major_axis_au: f64,
) -> (f64, f64) {
    let true_anomaly = 2.0
        * ((1.0 + eccentricity).sqrt() * (eccentric_anomaly * 0.5).tan()
            / (1.0 - eccentricity).sqrt())
        .atan();
    let radius = semi_major_axis_au * (1.0 - eccentricity * eccentric_anomaly.cos());
    (true_anomaly, radius)
}

/// Calculates a heliocentric position for an orbit whose explicit mean motion is
/// authoritative.
pub fn calculate_mean_motion_position(
    orbit: &MeanMotionOrbit,
    julian_date: JulianDate,
) -> Result<EclipticMeanJ2000<AstronomicalUnit>, ConicError> {
    let geometry = orbit.validated_geometry()?;
    let semi_major_axis = geometry.shape.semi_major_axis.value();
    let eccentricity = geometry.shape.eccentricity;
    let orientation = geometry.orientation;

    let dt_days = (julian_date - orbit.epoch).value();
    let mean_anomaly_deg = orbit.mean_motion_deg_per_day * dt_days;
    let mean_anomaly = Degrees::new(mean_anomaly_deg).to::<Radian>();
    let eccentric_anomaly = solve_keplers_equation(mean_anomaly, eccentricity);
    let (true_anomaly, radius) =
        elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, semi_major_axis);

    Ok(rotate_to_ecliptic(
        radius,
        orientation.inclination,
        orientation.argument_of_periapsis,
        orientation.longitude_of_ascending_node,
        true_anomaly,
    ))
}

/// Calculates a heliocentric position for unified conic elements.
pub fn calculate_conic_position(
    orbit: &ConicOrbit,
    julian_date: JulianDate,
) -> Result<EclipticMeanJ2000<AstronomicalUnit>, ConicError> {
    let geometry = orbit.validated_geometry()?;
    let kind = geometry.kind().map_err(map_validation_error)?;
    let periapsis_distance = geometry.shape.periapsis_distance.value();
    let eccentricity = geometry.shape.eccentricity;
    let orientation = geometry.orientation;
    match kind {
        ConicKind::Elliptic => {
            let semi_major_axis = periapsis_distance / (1.0 - eccentricity);
            let mean_motion =
                GAUSSIAN_GRAVITATIONAL_CONSTANT / (semi_major_axis * semi_major_axis.sqrt());
            let dt_days = (julian_date - orbit.epoch).value();
            let mean_anomaly =
                orbit.mean_anomaly_at_epoch.to::<Radian>() + Radians::new(mean_motion * dt_days);
            let eccentric_anomaly = solve_keplers_equation(mean_anomaly, eccentricity);
            let (true_anomaly, radius) =
                elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, semi_major_axis);

            Ok(rotate_to_ecliptic(
                radius,
                orientation.inclination,
                orientation.argument_of_periapsis,
                orientation.longitude_of_ascending_node,
                true_anomaly,
            ))
        }
        ConicKind::Hyperbolic => {
            let semi_major_axis = periapsis_distance / (eccentricity - 1.0);
            let mean_motion =
                GAUSSIAN_GRAVITATIONAL_CONSTANT / (semi_major_axis * semi_major_axis.sqrt());
            let dt_days = (julian_date - orbit.epoch).value();
            let mean_anomaly =
                orbit.mean_anomaly_at_epoch.to::<Radian>().value() + mean_motion * dt_days;
            let hyperbolic_anomaly = solve_hyperbolic_anomaly(mean_anomaly, eccentricity);
            let true_anomaly = 2.0
                * ((eccentricity + 1.0).sqrt() * (hyperbolic_anomaly * 0.5).sinh())
                    .atan2((eccentricity - 1.0).sqrt() * (hyperbolic_anomaly * 0.5).cosh());
            let radius = semi_major_axis * (eccentricity * hyperbolic_anomaly.cosh() - 1.0);

            Ok(rotate_to_ecliptic(
                radius,
                orientation.inclination,
                orientation.argument_of_periapsis,
                orientation.longitude_of_ascending_node,
                true_anomaly,
            ))
        }
        ConicKind::Parabolic => Err(ConicError::ParabolicUnsupported),
    }
}

impl MeanMotionOrbit {
    /// Calculates heliocentric coordinates at a given Julian date.
    pub fn position_at(
        &self,
        julian_date: JulianDate,
    ) -> Result<EclipticMeanJ2000<AstronomicalUnit>, ConicError> {
        calculate_mean_motion_position(self, julian_date)
    }
}

impl ConicOrbit {
    /// Calculates heliocentric coordinates at a given Julian date.
    pub fn position_at(
        &self,
        julian_date: JulianDate,
    ) -> Result<EclipticMeanJ2000<AstronomicalUnit>, ConicError> {
        calculate_conic_position(self, julian_date)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::macros::assert_cartesian_eq;

    #[test]
    fn mean_motion_position_is_at_periapsis_at_epoch() {
        let orbit = MeanMotionOrbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            0.9856076686,
            JulianDate::J2000,
        );
        let position = orbit.position_at(JulianDate::J2000).unwrap();
        assert_cartesian_eq!(position, EclipticMeanJ2000::new(1.0, 0.0, 0.0), 1e-10);
    }

    #[test]
    fn hyperbolic_position_is_finite() {
        let orbit = ConicOrbit::new(
            0.43355636 * AU,
            1.000956094769503,
            Degrees::new(110.804751),
            Degrees::new(260.04078),
            Degrees::new(68.15913068),
            Degrees::new(0.0),
            JulianDate::new(2_458_997.030_358_636_3),
        );
        let position = orbit.position_at(JulianDate::new(2458999.0)).unwrap();
        assert!(position.x().value().is_finite());
        assert!(position.y().value().is_finite());
        assert!(position.z().value().is_finite());
    }

    #[test]
    fn parabolic_orbit_is_rejected() {
        let orbit = ConicOrbit::new(
            1.0 * AU,
            1.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        assert!(matches!(
            orbit.position_at(JulianDate::J2000),
            Err(ConicError::ParabolicUnsupported)
        ));
    }

    #[test]
    fn invalid_mean_motion_semi_major_axis_maps_to_existing_error() {
        let orbit = MeanMotionOrbit::new(
            -1.0 * AU,
            0.1,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            1.0,
            JulianDate::J2000,
        );
        assert!(matches!(
            orbit.position_at(JulianDate::J2000),
            Err(ConicError::InvalidSemiMajorAxis)
        ));
    }

    #[test]
    fn invalid_mean_motion_eccentricity_maps_to_existing_error() {
        let orbit = MeanMotionOrbit::new(
            1.0 * AU,
            1.1,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            1.0,
            JulianDate::J2000,
        );
        assert!(matches!(
            orbit.position_at(JulianDate::J2000),
            Err(ConicError::InvalidEccentricity)
        ));
    }

    #[test]
    fn invalid_periapsis_distance_maps_to_existing_error() {
        let orbit = ConicOrbit::new(
            0.0 * AU,
            0.5,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        assert!(matches!(
            orbit.position_at(JulianDate::J2000),
            Err(ConicError::InvalidPeriapsisDistance)
        ));
    }
}
