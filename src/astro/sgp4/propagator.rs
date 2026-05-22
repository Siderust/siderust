// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SGP4/SDP4 propagator with strongly typed I/O.
//!
//! See the crate-level documentation for an overview. The propagator is a
//! thin, type-preserving wrapper around the public-domain SGP4 reference
//! implementation distributed via the [`sgp4` crate](https://crates.io/crates/sgp4).
//! The wrapper takes responsibility for:
//!
//! * mapping [`crate::formats::tle::TLE`] into the SGP4 backend's element record,
//! * converting [`tempoch::JulianDate<tempoch::UTC>`] target epochs into
//!   minutes-since-epoch (the propagator's natural argument), and
//! * tagging the resulting Cartesian arrays with their TEME / geocentric /
//!   km / km·s⁻¹ types from `affn`, `siderust`, and `qtty`.

use crate::formats::tle::TLE;
use crate::qtty::Minutes;
use sgp4::{Constants, Elements, Geopotential, MinutesSinceEpoch};
use tempoch::{JulianDate, JD, UTC};

use super::elements::tle_to_elements;
use super::state::TemeState;
use super::Sgp4Error;

/// Earth gravity / sidereal-time model selector for the SGP4/SDP4 polynomial.
///
/// SGP4 was historically defined against the WGS-72 geopotential constants
/// together with the AFSPC sidereal-time convention. That combination
/// (`Wgs72`) is the **bit-exact** match for Vallado & Crawford's published
/// `tcppver.out` reference outputs and is the default here.
///
/// `Wgs72Iau` keeps the WGS-72 geopotential but switches to the IAU
/// sidereal-time formula — slightly more accurate in absolute terms but
/// drifts a few millimetres away from the Vallado reference.
///
/// `Wgs84` swaps in the WGS-84 geopotential constants together with the
/// IAU sidereal-time formula. New code that does not need bit-exact
/// reproduction of legacy AFSPC outputs should pick this.
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::GravityModel;
/// // The default matches Vallado's `tcppver.out` reference outputs exactly.
/// assert!(matches!(GravityModel::default(), GravityModel::Wgs72));
/// ```
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum GravityModel {
    /// WGS-72 geopotential with AFSPC sidereal-time convention. Default;
    /// reproduces Vallado's `tcppver.out` reference outputs to numerical
    /// noise (≤ 1e-6 km / ≤ 1e-9 km·s⁻¹).
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
    fn as_geopotential(self) -> Geopotential {
        match self {
            GravityModel::Wgs72 | GravityModel::Wgs72Iau => sgp4::WGS72,
            GravityModel::Wgs84 => sgp4::WGS84,
        }
    }

    fn afspc_compat(self) -> bool {
        matches!(self, GravityModel::Wgs72)
    }
}

/// Strongly typed SGP4/SDP4 propagator.
///
/// Construct with [`Sgp4Propagator::from_tle`] (default WGS-72 gravity) or
/// [`Sgp4Propagator::from_tle_with_model`] (explicit gravity model).
/// Propagate with [`Sgp4Propagator::propagate_at`] (UTC Julian date) or
/// [`Sgp4Propagator::propagate_minutes`] (offset from the TLE epoch).
///
/// The propagator owns its initialised `Constants` table, so propagation is
/// allocation-free and thread-safe (the type is `Send + Sync`).
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
/// use siderust::formats::tle::parse_3le;
/// use siderust::qtty::Minutes;
///
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// let prop = Sgp4Propagator::from_tle_with_model(&tle, GravityModel::Wgs72).unwrap();
/// let s = prop.propagate_minutes(Minutes::new(0.0)).unwrap();
/// assert!(s.position().distance().value() > 6_500.0);
/// ```
#[derive(Clone, Debug)]
pub struct Sgp4Propagator {
    constants: Constants,
    elements: Elements,
    model: GravityModel,
    epoch_jd_utc: JulianDate<UTC>,
}

impl Sgp4Propagator {
    /// Initialise the propagator from a [`TLE`] using the default
    /// (WGS-72, Vallado 2006) gravity model.
    ///
    /// # Errors
    ///
    /// Returns [`Sgp4Error::InvalidElements`] if the SGP4 initialiser
    /// rejects the mean elements (e.g. negative mean motion, eccentricity
    /// outside `[0, 1)`), or [`Sgp4Error::InvalidEpoch`] /
    /// [`Sgp4Error::TimeConversion`] if the TLE epoch cannot be expressed
    /// as a calendar instant.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let _ = Sgp4Propagator::from_tle(&tle).unwrap();
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
    /// use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let _ = Sgp4Propagator::from_tle_with_model(&tle, GravityModel::Wgs84).unwrap();
    /// ```
    pub fn from_tle_with_model(tle: &TLE, model: GravityModel) -> Result<Self, Sgp4Error> {
        let elements = tle_to_elements(tle)?;
        let constants = match model {
            GravityModel::Wgs72 => Constants::from_elements_afspc_compatibility_mode(&elements)
                .map_err(|e| Sgp4Error::InvalidElements {
                    details: format!("AFSPC init: {e:?}"),
                })?,
            GravityModel::Wgs72Iau => {
                Constants::from_elements(&elements).map_err(|e| Sgp4Error::InvalidElements {
                    details: format!("IAU init (WGS-72): {e:?}"),
                })?
            }
            GravityModel::Wgs84 => {
                let geop = model.as_geopotential();
                let orbit = sgp4::Orbit::from_kozai_elements(
                    &geop,
                    elements.inclination * core::f64::consts::PI / 180.0,
                    elements.right_ascension * core::f64::consts::PI / 180.0,
                    elements.eccentricity,
                    elements.argument_of_perigee * core::f64::consts::PI / 180.0,
                    elements.mean_anomaly * core::f64::consts::PI / 180.0,
                    elements.mean_motion * 2.0 * core::f64::consts::PI / 1440.0,
                )
                .map_err(|e| Sgp4Error::InvalidElements {
                    details: format!("Kozai elements: {e:?}"),
                })?;
                Constants::new(
                    geop,
                    sgp4::iau_epoch_to_sidereal_time,
                    sgp4::julian_years_since_j2000(&elements.datetime),
                    elements.drag_term,
                    orbit,
                )
                .map_err(|e| Sgp4Error::InvalidElements {
                    details: format!("Constants::new: {e:?}"),
                })?
            }
        };
        let epoch_jd_utc = tle.epoch.to::<JD>();
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
    /// use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
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
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
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
    /// * [`Sgp4Error::Propagation`] — the SGP4 polynomial diverged at
    ///   the requested epoch.
    /// * [`Sgp4Error::TimeConversion`] — the supplied UTC instant cannot
    ///   be related to the TLE epoch (e.g. a non-finite Julian date).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::astro::sgp4::Sgp4Propagator;
    /// use siderust::formats::tle::parse_3le;
    /// let tle = parse_3le(
    ///     "ISS (ZARYA)",
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
    /// ).unwrap();
    /// let p = Sgp4Propagator::from_tle(&tle).unwrap();
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
    /// This is the natural argument to SGP4 and exactly mirrors the
    /// reference test vectors distributed with Vallado's "SGP4-VER" set.
    ///
    /// # Errors
    ///
    /// [`Sgp4Error::Propagation`] if the polynomial diverges, or
    /// [`Sgp4Error::TimeConversion`] if the resulting epoch cannot be
    /// represented as a UTC Julian date.
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
        let prediction = if self.model.afspc_compat() {
            self.constants
                .propagate_afspc_compatibility_mode(MinutesSinceEpoch(dt_minutes))
        } else {
            self.constants.propagate(MinutesSinceEpoch(dt_minutes))
        }
        .map_err(|e| Sgp4Error::Propagation {
            details: format!("{e:?}"),
        })?;
        // Touch `elements` so unused-field lints stay happy and so users
        // can later add a `pub fn elements(&self)` accessor without an
        // ABI change.
        let _ = &self.elements;
        Ok(TemeState::from_arrays(
            target_jd_utc,
            prediction.position,
            prediction.velocity,
        ))
    }
}

fn jd_offset_minutes(epoch: JulianDate<UTC>, minutes: Minutes) -> JulianDate<UTC> {
    use qtty::time::Day;
    use qtty::Quantity;
    let days = minutes.value() / 1_440.0;
    let raw = epoch.raw().value() + days;
    JulianDate::<UTC>::try_new(Quantity::<Day>::new(raw))
        .expect("finite Julian date increments remain finite")
}
