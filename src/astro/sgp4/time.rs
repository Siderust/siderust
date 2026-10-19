// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Nominal UTC calendar arithmetic for SGP4, independent of SI leap durations.
use super::Sgp4Error;
use crate::qtty::Minutes;
use alloc::string::ToString;
use chrono::{NaiveDateTime, Timelike};
use core::f64::consts::{PI, TAU};
#[cfg(not(feature = "std"))]
use qtty::{Real, Scalar, Transcendental};

pub(super) fn minutes_between(
    epoch: NaiveDateTime,
    target: NaiveDateTime,
) -> Result<Minutes, Sgp4Error> {
    target
        .signed_duration_since(epoch)
        .num_nanoseconds()
        .map(|ns| Minutes::new(ns as f64 / 60e9))
        .ok_or_else(|| {
            Sgp4Error::TimeConversion("nominal UTC offset exceeds nanosecond calendar range".into())
        })
}

pub(super) fn offset_epoch(
    epoch: NaiveDateTime,
    minutes: Minutes,
) -> Result<NaiveDateTime, Sgp4Error> {
    // Tempoch is the canonical finite/range validator for typed durations.
    let exact = tempoch::ExactDuration::try_from_quantity(minutes)
        .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?;
    let nanos = exact
        .as_nanoseconds_i()
        .map_err(|e| Sgp4Error::TimeConversion(e.to_string()))?
        .value();
    epoch
        .checked_add_signed(chrono::TimeDelta::nanoseconds(nanos))
        .ok_or_else(|| Sgp4Error::TimeConversion("UTC calendar overflow".into()))
}

/// Days since 1950 Jan 0 and epoch Greenwich sidereal angle. AFSPC preserves
/// the legacy scalar-JD rounding during initialization; IAU uses the split
/// calendar fraction. Neither expression is used for propagation offsets.
pub(super) fn model_epoch(datetime: NaiveDateTime, afspc: bool) -> (f64, f64) {
    let whole_seconds = datetime.and_utc().timestamp() - 946728000;
    let days = whole_seconds as f64 / 86400.0 + datetime.nanosecond() as f64 / 86400e9;
    let years = if afspc {
        ((2451545.0 + days) - 2451545.0) / 365.25
    } else {
        days / 365.25
    };
    let days_since_1950 = years * 365.25 + 18263.5;
    let sidereal = if afspc {
        let d = (years + 30.0) * 365.25 + 1.0;
        let integer = (d + 1e-8).floor();
        let fraction = d - integer;
        1.7321343856509374
            + 0.017202791694070364 * integer
            + (0.017202791694070364 + TAU) * fraction
            + d * d * 5.075_514_194_322_695e-15
    } else {
        let century = years / 100.0;
        (-6.2e-6 * century.powi(3)
            + 0.093104 * century * century
            + (876600.0 * 3600.0 + 8640184.812866) * century
            + 67310.54841)
            * (PI / 180.0)
            / 240.0
    };
    (days_since_1950, sidereal.rem_euclid(TAU))
}
