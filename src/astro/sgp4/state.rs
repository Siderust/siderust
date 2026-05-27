// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Strongly typed SGP4/SDP4 prediction state.
//!
//! [`TemeState`] bundles the position and velocity returned by a single
//! [`Sgp4Propagator::propagate_at`](crate::astro::sgp4::Sgp4Propagator::propagate_at) call,
//! together with the UTC instant the prediction is valid for.
//!
//! The state lives in the **TEME** (True Equator, Mean Equinox) frame defined
//! by `affn::frames::TEME`, geocentric, with positions in
//! [`qtty::length::Kilometer`] and velocities in [`KilometerPerSecond`].
//! Downstream conversion to ITRF / GCRF / ECEF is provided by
//! `siderust::coordinates::transform::providers::frames_teme` once Earth
//! Orientation Parameters are supplied.

use crate::coordinates::cartesian::{position, velocity};
use qtty::length::Kilometer;
use qtty::time::Second;
use qtty::Per;
use tempoch::{JulianDate, UTC};

/// Velocity unit alias used by the SGP4 propagator: kilometres per second.
///
/// SGP4 publishes velocities in km·s⁻¹; we expose them as a typed
/// [`qtty::Per<Kilometer, Second>`] so users cannot accidentally
/// mix them with m·s⁻¹ values from other parts of the stack.
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::KilometerPerSecond;
/// use qtty::Quantity;
/// let v: Quantity<KilometerPerSecond> = Quantity::new(7.5);
/// assert!((v.value() - 7.5).abs() < 1e-12);
/// ```
pub type KilometerPerSecond = Per<Kilometer, Second>;

/// Geocentric **TEME** Cartesian position with kilometre units.
///
/// Re-export of [`crate::coordinates::cartesian::position::TEME`] specialised
/// to [`qtty::length::Kilometer`].
pub type TemePositionKm = position::TEME<Kilometer>;

/// **TEME** Cartesian velocity with km·s⁻¹ units.
pub type TemeVelocityKmPerSec = velocity::TEME<KilometerPerSecond>;

/// SGP4 prediction at a single epoch in the **TEME** frame.
///
/// Construction is performed by the propagator; this struct only exposes
/// typed accessors. The fields are public-but-stable — adding new components
/// (e.g. covariance) would be a breaking change requiring a new type.
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::{Sgp4Propagator, TemeState};
/// use siderust::formats::tle::parse_3le;
/// use siderust::qtty::Minutes;
///
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// let prop = Sgp4Propagator::from_tle(&tle).unwrap();
/// let state: TemeState = prop.propagate_minutes(Minutes::new(0.0)).unwrap();
/// // Magnitudes are in kilometres and km/s for an LEO.
/// let r = state.position().as_array();
/// assert!(r[0].value().abs() < 8_000.0);
/// ```
#[derive(Clone, Debug)]
pub struct TemeState {
    epoch_utc: JulianDate<UTC>,
    position: TemePositionKm,
    velocity: TemeVelocityKmPerSec,
}

impl TemeState {
    /// Construct a typed state from the SGP4 raw arrays plus its valid epoch.
    ///
    /// `position_km` / `velocity_km_per_s` follow the standard SGP4 ordering
    /// `[x, y, z]` in TEME.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::TemeState;
    /// use tempoch::{JulianDate, Time, UTC};
    /// use chrono::{TimeZone, Utc};
    ///
    /// let t = Time::<UTC>::try_from_chrono(
    ///     Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap(),
    /// )
    /// .unwrap();
    /// let jd = t.to::<tempoch::JD>();
    /// let s = TemeState::from_arrays(
    ///     jd,
    ///     [7000.0, 0.0, 0.0],
    ///     [0.0, 7.5, 0.0],
    /// );
    /// assert!((s.position().x().value() - 7000.0).abs() < 1e-12);
    /// assert!((s.velocity().y().value() - 7.5).abs() < 1e-12);
    /// ```
    pub fn from_arrays(
        epoch_utc: JulianDate<UTC>,
        position_km: [f64; 3],
        velocity_km_per_s: [f64; 3],
    ) -> Self {
        use qtty::Quantity;
        let p = TemePositionKm::new(
            Quantity::<Kilometer>::new(position_km[0]),
            Quantity::<Kilometer>::new(position_km[1]),
            Quantity::<Kilometer>::new(position_km[2]),
        );
        let v = TemeVelocityKmPerSec::new(
            Quantity::<KilometerPerSecond>::new(velocity_km_per_s[0]),
            Quantity::<KilometerPerSecond>::new(velocity_km_per_s[1]),
            Quantity::<KilometerPerSecond>::new(velocity_km_per_s[2]),
        );
        Self {
            epoch_utc,
            position: p,
            velocity: v,
        }
    }

    /// UTC epoch the prediction is valid for.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// use siderust::qtty::Minutes;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
    /// let s = p.propagate_minutes(Minutes::new(0.0)).unwrap();
    /// assert!(s.epoch().raw().value() > 2_454_000.0);
    /// ```
    pub fn epoch(&self) -> JulianDate<UTC> {
        self.epoch_utc
    }

    /// Geocentric TEME position in kilometres.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// use siderust::qtty::Minutes;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
    /// let s = p.propagate_minutes(Minutes::new(0.0)).unwrap();
    /// let r = s.position();
    /// assert!(r.distance().value() > 6_000.0);
    /// ```
    pub fn position(&self) -> &TemePositionKm {
        &self.position
    }

    /// TEME velocity in km·s⁻¹.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// use siderust::qtty::Minutes;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
    /// let s = p.propagate_minutes(Minutes::new(0.0)).unwrap();
    /// let v = s.velocity();
    /// // LEO speed is ~7.7 km/s.
    /// assert!(v.magnitude().value() > 5.0 && v.magnitude().value() < 10.0);
    /// ```
    pub fn velocity(&self) -> &TemeVelocityKmPerSec {
        &self.velocity
    }
}
