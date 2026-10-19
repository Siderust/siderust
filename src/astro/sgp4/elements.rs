// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Adapt typed TLE mean elements to the native normalized numerical model.
//!
//! Callers always interact with the typed [`crate::formats::tle::TLE`]
//! record. The numerical model uses normalized Earth-radius/minute units.

use crate::formats::tle::TLE;
use alloc::format;
use chrono::Timelike;
use qtty::angular_rate::AngularRate;
use qtty::unit::{Minute, Radian};

use super::Sgp4Error;

/// Compact model inputs; angles are radians and mean motion is rad/min.
#[derive(Clone, Debug)]
pub(super) struct NativeElements {
    pub datetime: chrono::NaiveDateTime,
    pub orbit: super::model::Orbit,
    pub bstar: f64,
}

/// Convert a typed [`TLE`] into native propagation elements.
pub(crate) fn tle_to_elements(
    tle: &TLE,
    ctx: &tempoch::TimeContext,
) -> Result<NativeElements, Sgp4Error> {
    let dt_utc = tle
        .epoch
        .try_to_chrono_with(ctx)
        .map_err(|e| Sgp4Error::TimeConversion(format!("TLE epoch → chrono failed: {e:?}")))?;
    // SGP4 wants a NaiveDateTime in UTC.
    let datetime = dt_utc.naive_utc();
    if datetime.nanosecond() >= 1_000_000_000 {
        return Err(Sgp4Error::InvalidEpoch(
            "leap-second label is outside the SGP4 time grid",
        ));
    }
    if keplerian::Eccentricity::new_elliptic(tle.eccentricity).is_err() {
        return Err(Sgp4Error::InvalidElements {
            details: format!("eccentricity outside [0, 1): {}", tle.eccentricity),
        });
    }
    if !tle.mean_motion.value().is_finite() || tle.mean_motion.value() <= 0.0 {
        return Err(Sgp4Error::InvalidElements {
            details: format!("non-positive mean motion: {}", tle.mean_motion.value()),
        });
    }

    // Validate all independently mutable numerical inputs before initialization.
    if [
        tle.inclination.value(),
        tle.raan.value(),
        tle.argument_of_perigee.value(),
        tle.mean_anomaly.value(),
        tle.bstar,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(Sgp4Error::InvalidElements {
            details: "non-finite angle or BSTAR".into(),
        });
    }
    if !(0.0..=180.0).contains(&tle.inclination.value()) {
        return Err(Sgp4Error::InvalidElements {
            details: "inclination outside [0, 180] degrees".into(),
        });
    }
    // Preserve the legacy initializer's rejection of extreme mean motions:
    // the original multiply-by-PI normalization overflowed above this bound.
    if tle.mean_motion.value() > f64::MAX / core::f64::consts::PI {
        return Err(Sgp4Error::InvalidElements {
            details: "mean motion exceeds finite normalized range".into(),
        });
    }
    let mean_motion: AngularRate<Radian, Minute> = tle.mean_motion.to();
    Ok(NativeElements {
        datetime,
        orbit: super::model::Orbit {
            e: tle.eccentricity,
            i: tle.inclination.to::<qtty::unit::Radian>().value(),
            node: tle.raan.to::<qtty::unit::Radian>().value(),
            argp: tle.argument_of_perigee.to::<qtty::unit::Radian>().value(),
            m: tle.mean_anomaly.to::<qtty::unit::Radian>().value(),
            n: mean_motion.value(),
        },
        bstar: tle.bstar,
    })
}
