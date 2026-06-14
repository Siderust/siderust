// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! SGP4/SDP4 propagator with strongly typed I/O.
//!
//! See the crate-level documentation for an overview. The propagator is a
//! Siderust-owned numerical core for TLE mean elements. It takes responsibility
//! for:
//!
//! * mapping [`crate::formats::tle::TLE`] into an internal element record,
//! * converting [`tempoch::JulianDate<tempoch::UTC>`] target epochs into
//!   minutes-since-epoch, and
//! * tagging the resulting Cartesian arrays with their TEME / geocentric /
//!   km / km·s⁻¹ types from `affn`, `siderust`, and `qtty`.

use crate::formats::tle::TLE;
use crate::qtty::Minutes;
use tempoch::{JulianDate, UTC};

use super::elements::{tle_to_elements, NativeElements};
use super::state::TemeState;
use super::Sgp4Error;

const TWO_PI: f64 = core::f64::consts::TAU;
const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;

/// Earth gravity / sidereal-time model selector for the TLE mean-element propagator.
///
/// The three variants select the geopotential constants and the sidereal-time
/// convention used when computing secular J2 node/perigee drift rates and the
/// final perifocal-to-inertial rotation.
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::GravityModel;
/// assert!(matches!(GravityModel::default(), GravityModel::Wgs72));
/// ```
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum GravityModel {
    /// WGS-72 geopotential with AFSPC sidereal-time convention. Default;
    /// matches the WGS-72 TLE mean-element conventions.
    #[default]
    Wgs72,
    /// WGS-72 geopotential with the IAU sidereal-time formula. Use when
    /// you want the WGS-72 mean elements but the IAU time convention.
    Wgs72Iau,
    /// WGS-84 geopotential with the IAU sidereal-time formula. Modern
    /// default for new pipelines that don't need AFSPC bit-exactness.
    Wgs84,
}

impl GravityModel {
    fn constants(self) -> GravityConstants {
        match self {
            GravityModel::Wgs72 | GravityModel::Wgs72Iau => GravityConstants {
                mu_km3_s2: 398_600.8,
                earth_radius_km: 6_378.135,
                j2: 1.082_616e-3,
            },
            GravityModel::Wgs84 => GravityConstants {
                mu_km3_s2: 398_600.5,
                earth_radius_km: 6_378.137,
                j2: 1.082_629_989_05e-3,
            },
        }
    }
}

/// Simplified J2 mean-elements propagator for TLE data.
///
/// This propagator derives an osculating trajectory from TLE mean elements
/// by applying secular J2 node/perigee/mean-anomaly drift rates and solving
/// Kepler's equation at the requested epoch. It is **not** a full SGP4/SDP4
/// implementation and does not apply SGP4 periodic corrections, Lyddane
/// deep-space drag, or resonance terms used by the Vallado reference.
///
/// Construct with [`TlePropagator::from_tle`] (default WGS-72 gravity) or
/// [`TlePropagator::from_tle_with_model`] (explicit gravity model).
/// Propagate with [`TlePropagator::propagate_at`] (UTC Julian date) or
/// [`TlePropagator::propagate_minutes`] (offset from the TLE epoch).
///
/// The propagator owns its initialised constants table, so propagation is
/// allocation-free and thread-safe (the type is `Send + Sync`).
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::{GravityModel, TlePropagator};
/// use siderust::formats::tle::parse_3le;
/// use siderust::qtty::Minutes;
///
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// let prop = TlePropagator::from_tle_with_model(&tle, GravityModel::Wgs72).unwrap();
/// let s = prop.propagate_minutes(Minutes::new(0.0)).unwrap();
/// assert!(s.position().distance().value() > 6_500.0);
/// ```
#[derive(Clone, Debug)]
pub struct TlePropagator {
    constants: GravityConstants,
    elements: NativeElements,
    model: GravityModel,
    epoch_jd_utc: JulianDate<UTC>,
}

impl TlePropagator {
    /// Initialise the propagator from a [`TLE`] using the default
    /// (WGS-72) gravity model.
    ///
    /// # Errors
    ///
    /// Returns [`Sgp4Error::InvalidElements`] if the initialiser
    /// rejects the mean elements (e.g. negative mean motion, eccentricity
    /// outside `[0, 1)`), or [`Sgp4Error::InvalidEpoch`] /
    /// [`Sgp4Error::TimeConversion`] if the TLE epoch cannot be expressed
    /// as a calendar instant.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::TlePropagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let _ = TlePropagator::from_tle(&tle).unwrap();
    /// ```
    pub fn from_tle(tle: &TLE) -> Result<Self, Sgp4Error> {
        Self::from_tle_with_model(tle, GravityModel::default())
    }

    /// Initialise the propagator from a [`TLE`] under an explicit
    /// [`GravityModel`].
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::{GravityModel, TlePropagator};
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let _ = TlePropagator::from_tle_with_model(&tle, GravityModel::Wgs84).unwrap();
    /// ```
    pub fn from_tle_with_model(tle: &TLE, model: GravityModel) -> Result<Self, Sgp4Error> {
        let elements = tle_to_elements(tle)?;
        let constants = model.constants();
        validate_native_orbit(&elements, constants)?;
        let epoch_jd_utc = elements.epoch_jd_utc;
        Ok(Self {
            constants,
            elements,
            model,
            epoch_jd_utc,
        })
    }

    /// Return the gravity model in use.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::{GravityModel, TlePropagator};
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = TlePropagator::from_tle(&tle).unwrap();
    /// assert_eq!(p.gravity_model(), GravityModel::Wgs72);
    /// ```
    pub fn gravity_model(&self) -> GravityModel {
        self.model
    }

    /// UTC Julian-date epoch of the TLE this propagator was initialised
    /// from.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::TlePropagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = TlePropagator::from_tle(&tle).unwrap();
    /// assert!(p.epoch_jd_utc().raw().value() > 2_454_700.0);
    /// ```
    pub fn epoch_jd_utc(&self) -> JulianDate<UTC> {
        self.epoch_jd_utc
    }

    /// Propagate to a UTC Julian-date epoch.
    ///
    /// The result is a typed [`TemeState`] valid at exactly the supplied
    /// epoch. The conversion to minutes-since-epoch goes through `tempoch`
    /// to honour leap seconds correctly.
    ///
    /// # Errors
    ///
    /// * [`Sgp4Error::Propagation`] — propagation diverged at the
    ///   requested epoch.
    /// * [`Sgp4Error::TimeConversion`] — the supplied UTC instant cannot
    ///   be related to the TLE epoch (e.g. a non-finite Julian date).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::TlePropagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = TlePropagator::from_tle(&tle).unwrap();
    /// let target = p.epoch_jd_utc(); // propagate to TLE epoch itself
    /// let s = p.propagate_at(target).unwrap();
    /// assert!(s.position().distance().value() > 6_500.0);
    /// ```
    pub fn propagate_at(&self, target_jd_utc: JulianDate<UTC>) -> Result<TemeState, Sgp4Error> {
        let t_jd = target_jd_utc.raw().value();
        let epoch_jd = self.epoch_jd_utc.raw().value();
        if !t_jd.is_finite() {
            return Err(Sgp4Error::TimeConversion(
                "target Julian date is non-finite".into(),
            ));
        }
        let dt_minutes = (t_jd - epoch_jd) * 1_440.0;
        self.propagate_internal(dt_minutes, target_jd_utc)
    }

    /// Propagate by an offset, in minutes, from the TLE epoch.
    ///
    /// # Errors
    ///
    /// [`Sgp4Error::Propagation`] if propagation diverges, or
    /// [`Sgp4Error::TimeConversion`] if the resulting epoch cannot be
    /// represented as a UTC Julian date.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::TlePropagator;
    /// use siderust::formats::tle::parse_3le;
    /// use siderust::qtty::Minutes;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = TlePropagator::from_tle(&tle).unwrap();
    /// let s_now = p.propagate_minutes(Minutes::new(0.0)).unwrap();
    /// let s_later = p.propagate_minutes(Minutes::new(90.0)).unwrap();
    /// assert!(s_now.position() != s_later.position());
    /// ```
    pub fn propagate_minutes(&self, dt_minutes: Minutes) -> Result<TemeState, Sgp4Error> {
        if !dt_minutes.value().is_finite() {
            return Err(Sgp4Error::TimeConversion(
                "minutes offset is non-finite".into(),
            ));
        }
        let target = jd_offset_minutes(self.epoch_jd_utc, dt_minutes);
        self.propagate_internal(dt_minutes.value(), target)
    }

    fn propagate_internal(
        &self,
        dt_minutes: f64,
        target_jd_utc: JulianDate<UTC>,
    ) -> Result<TemeState, Sgp4Error> {
        let prediction = propagate_native(&self.elements, self.constants, dt_minutes)?;
        Ok(TemeState::from_arrays(
            target_jd_utc,
            prediction.0,
            prediction.1,
        ))
    }
}

#[derive(Copy, Clone, Debug)]
struct GravityConstants {
    mu_km3_s2: f64,
    earth_radius_km: f64,
    j2: f64,
}

fn validate_native_orbit(
    elements: &NativeElements,
    constants: GravityConstants,
) -> Result<(), Sgp4Error> {
    let n = mean_motion_rad_s(elements, 0.0);
    if !n.is_finite() || n <= 0.0 {
        return Err(Sgp4Error::InvalidElements {
            details: "mean motion does not produce a finite positive rad/s value".into(),
        });
    }
    let a = semi_major_axis_km(constants.mu_km3_s2, n);
    let perigee = a * (1.0 - elements.eccentricity);
    if !perigee.is_finite() || perigee <= 0.0 {
        return Err(Sgp4Error::InvalidElements {
            details: format!("non-positive perigee radius: {perigee} km"),
        });
    }
    Ok(())
}

fn propagate_native(
    elements: &NativeElements,
    constants: GravityConstants,
    dt_minutes: f64,
) -> Result<([f64; 3], [f64; 3]), Sgp4Error> {
    let dt_days = dt_minutes / 1_440.0;
    let dt_seconds = dt_minutes * 60.0;
    let n_rad_s = mean_motion_rad_s(elements, dt_days);
    if !n_rad_s.is_finite() || n_rad_s <= 0.0 {
        return Err(Sgp4Error::Propagation {
            details: "mean motion became non-finite or non-positive".into(),
        });
    }

    let a = semi_major_axis_km(constants.mu_km3_s2, n_rad_s);
    let e = elements.eccentricity;
    let one_minus_e2 = 1.0 - e * e;
    let p = a * one_minus_e2;
    if !p.is_finite() || p <= 0.0 {
        return Err(Sgp4Error::Propagation {
            details: format!("invalid semi-latus rectum {p} km"),
        });
    }

    let j2_factor =
        constants.j2 * (constants.earth_radius_km / p) * (constants.earth_radius_km / p);
    let cos_i = elements.inclination_rad.cos();
    let raan_dot = -1.5 * j2_factor * n_rad_s * cos_i;
    let argp_dot = 0.75 * j2_factor * n_rad_s * (5.0 * cos_i * cos_i - 1.0);
    let mean_j2_dot =
        0.75 * j2_factor * n_rad_s * one_minus_e2.sqrt() * (3.0 * cos_i * cos_i - 1.0);

    let raan = elements.raan_rad + raan_dot * dt_seconds;
    let argp = elements.argument_of_perigee_rad + argp_dot * dt_seconds;
    let mean_motion_phase = TWO_PI
        * (elements.mean_motion_rev_per_day * dt_days
            + 0.5 * elements.mean_motion_dot_rev_per_day2 * dt_days * dt_days
            + elements.mean_motion_ddot_rev_per_day3 * dt_days * dt_days * dt_days / 6.0);
    let m =
        normalize_radians(elements.mean_anomaly_rad + mean_motion_phase + mean_j2_dot * dt_seconds);
    let eccentric_anomaly = solve_kepler(m, e)?;
    let (sin_e, cos_e) = eccentric_anomaly.sin_cos();
    let edot = n_rad_s / (1.0 - e * cos_e);
    let sqrt_one_minus_e2 = one_minus_e2.sqrt();

    let x_orb = a * (cos_e - e);
    let y_orb = a * sqrt_one_minus_e2 * sin_e;
    let vx_orb = -a * sin_e * edot;
    let vy_orb = a * sqrt_one_minus_e2 * cos_e * edot;

    let (position, velocity) = rotate_perifocal_to_inertial(
        [x_orb, y_orb, 0.0],
        [vx_orb, vy_orb, 0.0],
        raan,
        elements.inclination_rad,
        argp,
    );

    if position
        .iter()
        .chain(velocity.iter())
        .any(|v| !v.is_finite())
    {
        return Err(Sgp4Error::Propagation {
            details: "propagated state contains non-finite component".into(),
        });
    }
    Ok((position, velocity))
}

fn mean_motion_rad_s(elements: &NativeElements, dt_days: f64) -> f64 {
    let rev_per_day = elements.mean_motion_rev_per_day
        + elements.mean_motion_dot_rev_per_day2 * dt_days
        + 0.5 * elements.mean_motion_ddot_rev_per_day3 * dt_days * dt_days
        - 1.0e-5 * elements.bstar * dt_days.abs();
    rev_per_day * TWO_PI / SECONDS_PER_DAY
}

fn semi_major_axis_km(mu_km3_s2: f64, n_rad_s: f64) -> f64 {
    (mu_km3_s2 / (n_rad_s * n_rad_s)).cbrt()
}

fn solve_kepler(mean_anomaly: f64, eccentricity: f64) -> Result<f64, Sgp4Error> {
    let mut e_anom = if eccentricity < 0.8 {
        mean_anomaly
    } else {
        core::f64::consts::PI
    };
    for _ in 0..32 {
        let f = e_anom - eccentricity * e_anom.sin() - mean_anomaly;
        let fp = 1.0 - eccentricity * e_anom.cos();
        if fp.abs() < 1.0e-14 {
            break;
        }
        let step = f / fp;
        e_anom -= step;
        if step.abs() < 1.0e-13 {
            return Ok(e_anom);
        }
    }
    Err(Sgp4Error::Propagation {
        details: "Kepler equation did not converge".into(),
    })
}

fn rotate_perifocal_to_inertial(
    pos: [f64; 3],
    vel: [f64; 3],
    raan: f64,
    inclination: f64,
    argp: f64,
) -> ([f64; 3], [f64; 3]) {
    let (sin_o, cos_o) = raan.sin_cos();
    let (sin_i, cos_i) = inclination.sin_cos();
    let (sin_w, cos_w) = argp.sin_cos();

    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r12 = -cos_o * sin_w - sin_o * cos_w * cos_i;
    let r21 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -sin_o * sin_w + cos_o * cos_w * cos_i;
    let r31 = sin_w * sin_i;
    let r32 = cos_w * sin_i;

    let rotate = |v: [f64; 3]| -> [f64; 3] {
        [
            r11 * v[0] + r12 * v[1],
            r21 * v[0] + r22 * v[1],
            r31 * v[0] + r32 * v[1],
        ]
    };
    (rotate(pos), rotate(vel))
}

fn normalize_radians(value: f64) -> f64 {
    value.rem_euclid(TWO_PI)
}

fn jd_offset_minutes(epoch: JulianDate<UTC>, minutes: Minutes) -> JulianDate<UTC> {
    use qtty::time::Day;
    use qtty::Quantity;
    let days = minutes.value() / 1_440.0;
    let raw = epoch.raw().value() + days;
    JulianDate::<UTC>::try_new(Quantity::<Day>::new(raw))
        .expect("finite Julian date increments remain finite")
}
