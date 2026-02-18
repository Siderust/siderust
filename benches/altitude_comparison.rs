// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Comparative benchmarks across body types for altitude computations.
//!
//! Measures single-point altitude evaluation and period-finding performance
//! for Sun, Moon, and Star (Sirius) side-by-side so that Criterion's HTML
//! reports show the relative cost of each body's altitude engine.
//!
//! Run with:
//! ```bash
//! cargo bench --bench altitude_comparison
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use qtty::*;
use siderust::bodies::{Moon, Sun};
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::direction;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period, MJD};
use std::hint::black_box;
use std::time::Duration;

fn build_period(days: u32) -> Period<MJD> {
    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap())
        + chrono::Duration::days(days as i64);

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    Period::new(
        ModifiedJulianDate::from_utc(start_dt),
        ModifiedJulianDate::from_utc(end_dt),
    )
}

fn sirius_icrs() -> direction::ICRS {
    direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716))
}

// =============================================================================
// Single-Point Altitude Evaluation
// =============================================================================

fn bench_single_altitude(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
    let mjd = ModifiedJulianDate::new(51544.5); // J2000
    let sirius = sirius_icrs();

    let mut group = c.benchmark_group("altitude/single_eval");

    group.bench_function("sun", |b| {
        b.iter(|| Sun.altitude_at(black_box(&site), black_box(mjd)));
    });

    group.bench_function("moon", |b| {
        b.iter(|| Moon.altitude_at(black_box(&site), black_box(mjd)));
    });

    group.bench_function("star_sirius", |b| {
        b.iter(|| sirius.altitude_at(black_box(&site), black_box(mjd)));
    });

    group.finish();
}

// =============================================================================
// 7-Day Above-Horizon Search
// =============================================================================

fn bench_above_horizon_7day(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
    let period = build_period(7);
    let sirius = sirius_icrs();

    let mut group = c.benchmark_group("altitude/above_horizon_7day");

    group.bench_function("sun", |b| {
        b.iter(|| {
            Sun.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.bench_function("moon", |b| {
        b.iter(|| {
            Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.bench_function("star_sirius", |b| {
        b.iter(|| {
            sirius.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.finish();
}

// =============================================================================
// 30-Day Below-Threshold Search
// =============================================================================

fn bench_below_threshold_30day(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
    let period = build_period(30);
    let sirius = sirius_icrs();

    let mut group = c.benchmark_group("altitude/below_threshold_30day");

    group.bench_function("sun_astro_night", |b| {
        b.iter(|| {
            Sun.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(-18.0)),
            )
        });
    });

    group.bench_function("moon_below_horizon", |b| {
        b.iter(|| {
            Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.bench_function("star_sirius_below_horizon", |b| {
        b.iter(|| {
            sirius.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.finish();
}

// =============================================================================
// 365-Day Period Search (full year — main performance target)
// =============================================================================

fn bench_period_search_365day(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
    let period = build_period(365);
    let sirius = sirius_icrs();

    let mut group = c.benchmark_group("altitude/period_search_365day");

    group.bench_function("sun_above_horizon", |b| {
        b.iter(|| {
            Sun.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.bench_function("moon_above_horizon", |b| {
        b.iter(|| {
            Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.bench_function("star_sirius_above_horizon", |b| {
        b.iter(|| {
            sirius.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            )
        });
    });

    group.finish();
}

criterion_group! {
    name = altitude_comparison_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(10))
        .sample_size(20)
        .without_plots();
    targets = bench_single_altitude,
              bench_above_horizon_7day,
              bench_below_threshold_30day,
              bench_period_search_365day
}
criterion_main!(altitude_comparison_benches);
