// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! SGP4/SDP4 propagator with strongly typed I/O.
//!
//! See the crate-level documentation for an overview. The propagator is a
//! native Vallado numerical model. It takes responsibility
//! for:
//!
//! * mapping [`crate::formats::tle::TLE`] into an internal element record,
//! * converting [`tempoch::JulianDate<tempoch::UTC>`] target epochs into
//!   minutes-since-epoch, and
//! * tagging the resulting Cartesian arrays with their TEME / geocentric /
//!   km / km·s⁻¹ types from `affn`, `siderust`, and `qtty`.

use crate::formats::tle::TLE;
use crate::qtty::Minutes;
use alloc::string::ToString;
use tempoch::{JulianDate, TimeContext, UTC};

use super::elements::tle_to_elements;
use super::model::Constants;
use super::state::TemeState;
use super::Sgp4Error;

/// Earth gravity / sidereal-time model selector for the TLE mean-element propagator.
///
/// The three variants select the geopotential constants and the sidereal-time
/// convention used by the full SGP4/SDP4 model. WGS-72/AFSPC is the actual
/// default and matches Vallado verification vectors; WGS-72/IAU and WGS-84/IAU
/// are explicit alternatives and can produce different results.
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
    /// WGS-84 geopotential with the IAU sidereal-time formula. Explicit
    /// alternative; this is not the default.
    Wgs84,
}

/// Full Vallado/AFSPC SGP4/SDP4 propagation of TLE mean elements.
///
/// Includes drag, periodic corrections, lunar/solar perturbations and deep-space
/// resonances in the native Siderust numerical model.
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
    context: TimeContext,
    constants: Constants,
    epoch_calendar_utc: chrono::NaiveDateTime,
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

    /// Initialise using a caller-provided, owned Tempoch time-data snapshot
    /// and the default WGS-72 gravity model.
    pub fn from_tle_with_context(tle: &TLE, context: TimeContext) -> Result<Self, Sgp4Error> {
        Self::from_tle_with_model_and_context(tle, GravityModel::default(), context)
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
        Self::from_tle_with_model_and_context(tle, model, TimeContext::new())
    }

    /// Initialise with an explicit gravity model and a caller-owned Tempoch
    /// time-data snapshot, retained for all subsequent epoch conversions.
    pub fn from_tle_with_model_and_context(
        tle: &TLE,
        model: GravityModel,
        context: TimeContext,
    ) -> Result<Self, Sgp4Error> {
        let elements = tle_to_elements(tle, &context)?;
        let constants = Constants::new(&elements, model)?;
        let epoch_jd_utc = tle
            .epoch
            .to_with::<tempoch::JD>(&context)
            .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?;
        Ok(Self {
            context,
            constants,
            epoch_calendar_utc: elements.datetime,
            model,
            epoch_jd_utc,
        })
    }

    /// Immutable time-data snapshot used for civil UTC encoding and decoding.
    pub fn time_context(&self) -> &TimeContext {
        &self.context
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
    /// assert!(p.epoch_jd_utc().try_raw_with(p.time_context()).unwrap().value() > 2_454_700.0);
    /// ```
    pub fn epoch_jd_utc(&self) -> JulianDate<UTC> {
        self.epoch_jd_utc
    }

    /// Propagate to a UTC Julian-date epoch.
    ///
    /// The result is a typed [`TemeState`] valid at exactly the supplied
    /// epoch. The model uses civil UTC days of 1,440 nominal minutes,
    /// not elapsed SI time. Leap-second instants cannot be encoded on this grid.
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
        // Validate the UTC grid first: leap-second labels have no nominal JD.
        target_jd_utc
            .try_raw_with(&self.context)
            .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?;
        let datetime = target_jd_utc
            .try_to_chrono_with(&self.context)
            .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?
            .naive_utc();
        let minutes = super::time::minutes_between(self.epoch_calendar_utc, datetime)?;
        self.propagate_internal(minutes, target_jd_utc)
    }

    /// Propagate by nominal model minutes (1,440 per civil UTC day) from the TLE epoch.
    ///
    /// # Errors
    ///
    /// [`Sgp4Error::Propagation`] if propagation diverges, or
    /// [`Sgp4Error::TimeConversion`] if the resulting epoch cannot be
    /// represented as a UTC Julian date. Offsets are limited to approximately
    /// ±292 years by the calendar conversion (nanosecond resolution).
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
        let datetime = super::time::offset_epoch(self.epoch_calendar_utc, dt_minutes)?;
        let target = tempoch::Time::<UTC>::try_from_chrono_with(datetime.and_utc(), &self.context)
            .and_then(|t| t.to_with::<tempoch::JD>(&self.context))
            .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?;
        self.propagate_internal(dt_minutes, target)
    }

    fn propagate_internal(
        &self,
        dt_minutes: Minutes,
        target_jd_utc: JulianDate<UTC>,
    ) -> Result<TemeState, Sgp4Error> {
        let (position, velocity) = self.constants.propagate(dt_minutes)?;
        Ok(TemeState::from_arrays(target_jd_utc, position, velocity))
    }
}
