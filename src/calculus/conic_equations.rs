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

use crate::astro::conic::{ConicError, ConicKind, ConicOrbit, MeanMotionOrbit};
use crate::astro::orbit::{KeplerianOrbit, OrientationTrig, PreparedOrbit};
use crate::astro::units::heliocentric_period_days;
use crate::coordinates::cartesian::position::EclipticMeanJ2000;
use crate::qtty::*;
use crate::time::JulianDate;
use keplerian::anomaly::{eccentric_from_mean, hyperbolic_from_mean, AnomalyOptions};
use keplerian::Eccentricity;

/// Gaussian gravitational constant `k` in AU^{3/2} d^{-1}.
///
/// Encodes `√(GM_☉)` in the heliocentric AU-day system. Using this
/// constant ties all propagation in this module to a Sun-centered
/// context.
const GAUSSIAN_GRAVITATIONAL_CONSTANT: f64 = 0.01720209895;

/// Converts a gravitational parameter from km³/s² to AU³/day².
///
/// The Gaussian constant satisfies k² = GM_☉ [AU³/day²], so:
///
/// ```text
/// mu [AU³/day²] = mu [km³/s²] × (86400 s/day)² / (149597870.7 km/AU)³
/// ```
#[inline]
fn mu_to_au3_d2(mu: GravitationalParameter) -> f64 {
    const KM_PER_AU: f64 = 149_597_870.7;
    const S_PER_DAY: f64 = 86_400.0;
    mu.value() * (S_PER_DAY * S_PER_DAY) / (KM_PER_AU * KM_PER_AU * KM_PER_AU)
}
/// Siderust keeps heliocentric orbit wrappers here, while `keplerian` owns the
/// reusable Kepler-equation solver.
#[inline]
fn solve_elliptic_anomaly(mean_anomaly: Radians, eccentricity: f64) -> Radians {
    let options = AnomalyOptions {
        max_iter: 100,
        tol: 1e-15,
    };
    eccentric_from_mean(
        mean_anomaly,
        Eccentricity::new_unchecked(eccentricity),
        options,
    )
    .expect("validated elliptic orbit must solve Kepler's equation")
}

fn solve_hyperbolic_anomaly(mean_anomaly_radians: f64, eccentricity: f64) -> Option<f64> {
    hyperbolic_from_mean(
        Radians::new(mean_anomaly_radians),
        Eccentricity::new_unchecked(eccentricity),
        AnomalyOptions {
            max_iter: 100,
            tol: 1e-14,
        },
    )
    .ok()
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
    let (sin_i, cos_i) = inclination.sin_cos();
    let (sin_u, cos_u) =
        (argument_of_periapsis.to::<Radian>() + Radians::new(true_anomaly_radians)).sin_cos();
    let (sin_node, cos_node) = longitude_of_ascending_node.sin_cos();
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
    let geometry = orbit.geometry();
    let semi_major_axis = geometry.shape().semi_major_axis().value();
    let eccentricity = geometry.shape().eccentricity();
    let orientation = geometry.orientation();

    let trig = OrientationTrig::from_orientation(orientation);
    let dt_days = (julian_date.raw() - orbit.epoch.raw()).value();
    let mean_anomaly_rad =
        (orbit.mean_motion.value().to_radians() * dt_days).rem_euclid(std::f64::consts::TAU);
    let mean_anomaly = Radians::new(mean_anomaly_rad);
    let eccentric_anomaly = solve_elliptic_anomaly(mean_anomaly, eccentricity);
    let (true_anomaly, radius) =
        elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, semi_major_axis);

    Ok(rotate_to_ecliptic_precomputed(radius, &trig, true_anomaly))
}

/// Calculates a heliocentric position for unified conic elements.
pub fn calculate_conic_position(
    orbit: &ConicOrbit,
    julian_date: JulianDate,
) -> Result<EclipticMeanJ2000<AstronomicalUnit>, ConicError> {
    let geometry = orbit.geometry();
    let kind = geometry.kind();
    let periapsis_distance = geometry.shape().periapsis_distance().value();
    let eccentricity = geometry.shape().eccentricity();
    let orientation = geometry.orientation();
    let trig = OrientationTrig::from_orientation(orientation);
    match kind {
        ConicKind::Elliptic => {
            let semi_major_axis = periapsis_distance / (1.0 - eccentricity);
            let mean_motion =
                GAUSSIAN_GRAVITATIONAL_CONSTANT / (semi_major_axis * semi_major_axis.sqrt());
            let dt_days = (julian_date.raw() - orbit.epoch.raw()).value();
            let mean_anomaly_raw =
                orbit.mean_anomaly_at_epoch.to::<Radian>().value() + mean_motion * dt_days;
            let mean_anomaly = Radians::new(mean_anomaly_raw.rem_euclid(std::f64::consts::TAU));
            let eccentric_anomaly = solve_elliptic_anomaly(mean_anomaly, eccentricity);
            let (true_anomaly, radius) =
                elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, semi_major_axis);

            Ok(rotate_to_ecliptic_precomputed(radius, &trig, true_anomaly))
        }
        ConicKind::Hyperbolic => {
            let semi_major_axis = periapsis_distance / (eccentricity - 1.0);
            let mean_motion =
                GAUSSIAN_GRAVITATIONAL_CONSTANT / (semi_major_axis * semi_major_axis.sqrt());
            let dt_days = (julian_date.raw() - orbit.epoch.raw()).value();
            let mean_anomaly =
                orbit.mean_anomaly_at_epoch.to::<Radian>().value() + mean_motion * dt_days;
            let hyperbolic_anomaly = solve_hyperbolic_anomaly(mean_anomaly, eccentricity)
                .ok_or(ConicError::HyperbolicSolverFailed)?;
            let true_anomaly = 2.0
                * ((eccentricity + 1.0).sqrt() * (hyperbolic_anomaly * 0.5).sinh())
                    .atan2((eccentricity - 1.0).sqrt() * (hyperbolic_anomaly * 0.5).cosh());
            let radius = semi_major_axis * (eccentricity * hyperbolic_anomaly.cosh() - 1.0);

            Ok(rotate_to_ecliptic_precomputed(radius, &trig, true_anomaly))
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

/// Calculates the heliocentric coordinates of a Keplerian orbit at a Julian date.
///
/// Uses the Gaussian gravitational constant (`k = 0.01720209895 AU^{3/2} d^{-1}`)
/// which encodes `GM_☉` in the AU-day system. Elements must be heliocentric and
/// in AU. For orbits around other central bodies use
/// [`calculate_orbit_position_with_mu`] with an explicit `mu`.
pub fn calculate_orbit_position(
    elements: &KeplerianOrbit,
    julian_date: JulianDate,
) -> EclipticMeanJ2000<AstronomicalUnit> {
    // Mean motion derived the same way as PreparedOrbit to ensure identical results.
    let a = elements.shape().semi_major_axis().value();
    type RadiansPerDay = crate::qtty::angular_rate::AngularRate<Radian, Day>;
    let period_days = heliocentric_period_days(a);
    let n = RadiansPerDay::new(std::f64::consts::TAU / period_days);
    let dt: Days = julian_date.raw() - elements.epoch.raw();
    let m0_rad = elements.mean_anomaly_at_epoch.to::<Radian>();
    let eccentricity = elements.shape().eccentricity();
    let mean_anomaly = (m0_rad + (n * dt).to::<Radian>()) % std::f64::consts::TAU;
    let eccentric_anomaly = solve_elliptic_anomaly(mean_anomaly, eccentricity);
    let (true_anomaly, radius) =
        elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, a);
    rotate_to_ecliptic(
        radius,
        elements.orientation().inclination(),
        elements.orientation().argument_of_periapsis(),
        elements.orientation().longitude_of_ascending_node(),
        true_anomaly,
    )
}

/// Calculates the position of a Keplerian orbit at a Julian date using an
/// explicit gravitational parameter.
///
/// The elements are assumed to be in AU, and the returned position is in the
/// ecliptic mean J2000 frame. `mu` may be the gravitational parameter of any
/// central body (not just the Sun): use [`GM_SUN`] for heliocentric orbits,
/// [`GM_EARTH`] for geocentric, or any other body's μ.
pub fn calculate_orbit_position_with_mu(
    elements: &KeplerianOrbit,
    julian_date: JulianDate,
    mu: GravitationalParameter,
) -> EclipticMeanJ2000<AstronomicalUnit> {
    let a_au = elements.shape().semi_major_axis().value();
    let mu_au3_d2 = mu_to_au3_d2(mu);
    type RadiansPerDay = crate::qtty::angular_rate::AngularRate<Radian, Day>;
    // mean motion n = sqrt(mu / a³) rad/day
    let n: RadiansPerDay = RadiansPerDay::new((mu_au3_d2 / (a_au * a_au * a_au)).sqrt());
    let dt: Days = julian_date.raw() - elements.epoch.raw();
    let m0_rad = elements.mean_anomaly_at_epoch.to::<Radian>();
    let eccentricity = elements.shape().eccentricity();
    let mean_anomaly = (m0_rad + (n * dt).to::<Radian>()) % std::f64::consts::TAU;
    let eccentric_anomaly = solve_elliptic_anomaly(mean_anomaly, eccentricity);
    let (true_anomaly, radius) =
        elliptic_true_anomaly_and_radius(eccentric_anomaly, eccentricity, a_au);

    rotate_to_ecliptic(
        radius,
        elements.orientation().inclination(),
        elements.orientation().argument_of_periapsis(),
        elements.orientation().longitude_of_ascending_node(),
        true_anomaly,
    )
}

impl KeplerianOrbit {
    /// Calculates heliocentric coordinates at a given Julian date.
    ///
    /// Uses the Gaussian gravitational constant, which encodes `GM_☉` in the
    /// AU-day system. This method is **heliocentric-only**. For orbits around
    /// any other central body supply the body's μ via [`position_with_mu`].
    ///
    /// [`position_with_mu`]: KeplerianOrbit::position_with_mu
    pub fn kepler_position(&self, jd: JulianDate) -> EclipticMeanJ2000<AstronomicalUnit> {
        calculate_prepared_position(&PreparedOrbit::from_validated(*self), jd)
    }

    /// Calculates orbital position at a given Julian date using an explicit
    /// gravitational parameter.
    ///
    /// Use this method when the orbit is not heliocentric. Supply the central
    /// body's μ from the `qtty` dynamics constants:
    ///
    /// - [`GM_SUN`] for heliocentric (same as [`kepler_position`])
    /// - [`GM_EARTH`] for geocentric (Earth-orbiting satellites, Moon)
    /// - [`GM_JUPITER`] for Jupiter-centric (Galilean moons)
    /// - etc.
    ///
    /// [`kepler_position`]: KeplerianOrbit::kepler_position
    pub fn position_with_mu(
        &self,
        jd: JulianDate,
        mu: GravitationalParameter,
    ) -> EclipticMeanJ2000<AstronomicalUnit> {
        calculate_orbit_position_with_mu(self, jd, mu)
    }
}

/// Fast propagation for a [`PreparedOrbit`] using its precomputed constants.
pub fn calculate_prepared_position(
    prepared: &PreparedOrbit,
    julian_date: JulianDate,
) -> EclipticMeanJ2000<AstronomicalUnit> {
    let dt = (julian_date.raw() - prepared.elements().epoch.raw()).value();
    let mean_anomaly = Radians::new(
        (prepared.m0_rad() + prepared.mean_motion().value() * dt) % std::f64::consts::TAU,
    );
    let eccentricity = prepared.elements().shape().eccentricity();
    let eccentric_anomaly = solve_elliptic_anomaly(mean_anomaly, eccentricity);
    let (true_anomaly, radius) = elliptic_true_anomaly_and_radius(
        eccentric_anomaly,
        eccentricity,
        prepared.elements().shape().semi_major_axis().value(),
    );

    rotate_to_ecliptic_precomputed(radius, prepared.orientation_trig(), true_anomaly)
}

impl PreparedOrbit {
    /// Calculates heliocentric coordinates at a given Julian date.
    #[inline]
    pub fn position_at(&self, jd: JulianDate) -> EclipticMeanJ2000<AstronomicalUnit> {
        calculate_prepared_position(self, jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::macros::assert_cartesian_eq;
    use crate::qtty::{angular_rate::AngularRate, Days};

    #[test]
    fn mean_motion_position_is_at_periapsis_at_epoch() {
        let orbit = MeanMotionOrbit::try_new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            AngularRate::<Degree, Day>::new(0.9856076686),
            crate::J2000,
        )
        .unwrap();
        let position = orbit.position_at(crate::J2000).unwrap();
        assert_cartesian_eq!(position, EclipticMeanJ2000::new(1.0, 0.0, 0.0), 1e-10);
    }

    #[test]
    fn mean_motion_large_dt_is_finite() {
        let orbit = MeanMotionOrbit::try_new(
            1.0 * AU,
            0.2,
            Degrees::new(10.0),
            Degrees::new(20.0),
            Degrees::new(30.0),
            AngularRate::<Degree, Day>::new(0.9856076686),
            crate::J2000,
        )
        .unwrap();
        // 1e8 days ~ 274,000 years — tests M normalization for large accumulation.
        let position = orbit
            .position_at(crate::time::JulianDate::new(2451545.0 + 1e8))
            .unwrap();
        assert!(position.x().is_finite());
        assert!(position.y().is_finite());
        assert!(position.z().is_finite());
    }

    #[test]
    fn hyperbolic_position_is_finite() {
        let orbit = ConicOrbit::try_new(
            0.43355636 * AU,
            1.000956094769503,
            Degrees::new(110.804751),
            Degrees::new(260.04078),
            Degrees::new(68.15913068),
            Degrees::new(0.0),
            crate::time::JulianDate::new(2_458_997.030_358_636_3),
        )
        .unwrap();
        let position = orbit
            .position_at(crate::time::JulianDate::new(2458999.0))
            .unwrap();
        assert!(position.x().is_finite());
        assert!(position.y().is_finite());
        assert!(position.z().is_finite());
    }

    #[test]
    fn parabolic_orbit_is_rejected_at_construction() {
        assert!(matches!(
            ConicOrbit::try_new(
                1.0 * AU,
                1.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                crate::J2000,
            ),
            Err(ConicError::ParabolicUnsupported)
        ));
    }

    #[test]
    fn invalid_mean_motion_semi_major_axis_maps_to_existing_error() {
        assert!(matches!(
            MeanMotionOrbit::try_new(
                -1.0 * AU,
                0.1,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(1.0),
                crate::J2000,
            ),
            Err(ConicError::InvalidSemiMajorAxis)
        ));
    }

    #[test]
    fn invalid_mean_motion_eccentricity_maps_to_hyperbolic_error() {
        assert!(matches!(
            MeanMotionOrbit::try_new(
                1.0 * AU,
                1.1,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(1.0),
                crate::J2000,
            ),
            Err(ConicError::HyperbolicNotSupported)
        ));
    }

    #[test]
    fn invalid_periapsis_distance_maps_to_existing_error() {
        assert!(matches!(
            ConicOrbit::try_new(
                0.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                crate::J2000,
            ),
            Err(ConicError::InvalidPeriapsisDistance)
        ));
    }

    #[test]
    fn near_parabolic_hyperbolic_position_is_finite() {
        // This is the reported failure case: e barely above 1 with small M.
        let orbit = ConicOrbit::try_new(
            1.0 * AU,
            1.0000001,
            Degrees::new(10.0),
            Degrees::new(20.0),
            Degrees::new(30.0),
            Degrees::new(0.001),
            crate::J2000,
        )
        .unwrap();
        let position = orbit
            .position_at(crate::time::JulianDate::new(2451545.5))
            .unwrap();
        assert!(position.x().is_finite());
        assert!(position.y().is_finite());
        assert!(position.z().is_finite());
    }

    #[test]
    fn moderate_hyperbolic_position_is_finite() {
        let orbit = ConicOrbit::try_new(
            1.0 * AU,
            2.0,
            Degrees::new(30.0),
            Degrees::new(60.0),
            Degrees::new(90.0),
            Degrees::new(0.0),
            crate::J2000,
        )
        .unwrap();
        let position = orbit
            .position_at(crate::time::JulianDate::new(2451645.0))
            .unwrap();
        assert!(position.x().is_finite());
        assert!(position.y().is_finite());
        assert!(position.z().is_finite());
    }

    #[test]
    fn large_mean_anomaly_hyperbolic_is_finite() {
        let orbit = ConicOrbit::try_new(
            1.0 * AU,
            1.5,
            Degrees::new(45.0),
            Degrees::new(90.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            crate::J2000,
        )
        .unwrap();
        // Large dt to produce large M
        let position = orbit
            .position_at(crate::time::JulianDate::new(2451545.0 + 10000.0))
            .unwrap();
        assert!(position.x().is_finite());
        assert!(position.y().is_finite());
        assert!(position.z().is_finite());
    }

    #[test]
    fn keplerian_orbit_position_is_at_periapsis_at_epoch() {
        let orbit = KeplerianOrbit::new(
            2.0 * AU,
            0.1,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            crate::J2000,
        );

        assert_cartesian_eq!(
            calculate_orbit_position(&orbit, crate::J2000),
            EclipticMeanJ2000::new(1.8, 0.0, 0.0),
            1e-10
        );
    }

    #[test]
    fn prepared_orbit_matches_keplerian_position() {
        let orbit = KeplerianOrbit::new(
            1.0 * AU,
            0.0167,
            Degrees::new(0.00005),
            Degrees::new(-11.26064),
            Degrees::new(102.94719),
            Degrees::new(100.46435),
            crate::J2000,
        );
        let jd = crate::time::JulianDate::new((crate::J2000.raw() + Days::new(100.0)).value());
        let prepared = PreparedOrbit::try_from(orbit).unwrap();

        let direct = orbit.kepler_position(jd);
        let cached = prepared.position_at(jd);
        assert!((direct.x() - cached.x()).abs() < AstronomicalUnits::new(1e-12));
        assert!((direct.y() - cached.y()).abs() < AstronomicalUnits::new(1e-12));
        assert!((direct.z() - cached.z()).abs() < AstronomicalUnits::new(1e-12));
    }
}
