// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Native TLE mean-element record used by the propagator.
//!
//! Callers always interact with the typed [`crate::formats::tle::TLE`]
//! record. This module copies the SGP4-relevant fields into a compact
//! crate-private representation with radians and revolutions-per-day scalars
//! for the numerical core.

use crate::formats::tle::TLE;
use chrono::Datelike;

use super::Sgp4Error;

/// Native copy of the TLE fields needed by the propagator.
///
/// Angles are stored in radians. Mean motion and derivatives keep the TLE
/// convention of revolutions per day and its time derivatives.
#[derive(Clone, Debug)]
pub(crate) struct NativeElements {
    pub(crate) epoch_jd_utc: tempoch::JulianDate<tempoch::UTC>,
    pub(crate) inclination_rad: f64,
    pub(crate) raan_rad: f64,
    pub(crate) eccentricity: f64,
    pub(crate) argument_of_perigee_rad: f64,
    pub(crate) mean_anomaly_rad: f64,
    pub(crate) mean_motion_rev_per_day: f64,
    pub(crate) mean_motion_dot_rev_per_day2: f64,
    pub(crate) mean_motion_ddot_rev_per_day3: f64,
    pub(crate) bstar: f64,
}

/// Convert a typed [`TLE`] into native propagation elements.
pub(crate) fn tle_to_elements(tle: &TLE) -> Result<NativeElements, Sgp4Error> {
    let dt_utc = tle
        .epoch
        .try_to_chrono()
        .map_err(|e| Sgp4Error::TimeConversion(format!("TLE epoch → chrono failed: {e:?}")))?;
    // SGP4 wants a NaiveDateTime in UTC.
    let datetime = dt_utc.naive_utc();
    // Sanity-check the year: AFSPC epoch encoding loses precision before 1957.
    if datetime.year() < 1950 {
        return Err(Sgp4Error::InvalidEpoch("year < 1950 not representable"));
    }
    if !tle.eccentricity.is_finite() || !(0.0..1.0).contains(&tle.eccentricity) {
        return Err(Sgp4Error::InvalidElements {
            details: format!("eccentricity outside [0, 1): {}", tle.eccentricity),
        });
    }
    if !tle.mean_motion.value().is_finite() || tle.mean_motion.value() <= 0.0 {
        return Err(Sgp4Error::InvalidElements {
            details: format!("non-positive mean motion: {}", tle.mean_motion.value()),
        });
    }

    Ok(NativeElements {
        epoch_jd_utc: tle.epoch.to::<tempoch::JD>(),
        inclination_rad: tle.inclination.value().to_radians(),
        raan_rad: tle.raan.value().to_radians(),
        eccentricity: tle.eccentricity,
        argument_of_perigee_rad: tle.argument_of_perigee.value().to_radians(),
        mean_anomaly_rad: tle.mean_anomaly.value().to_radians(),
        mean_motion_rev_per_day: tle.mean_motion.value(),
        mean_motion_dot_rev_per_day2: tle.mean_motion_dot,
        mean_motion_ddot_rev_per_day3: tle.mean_motion_ddot,
        bstar: tle.bstar,
    })
}
