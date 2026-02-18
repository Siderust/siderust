// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for lunar altitude calculations.
//!
//! Tests the performance of calculating Moon altitude and finding periods when
//! the Moon is above/below the horizon using the ELP2000 lunar theory.
//!
//! ## Algorithms Compared
//!
//! - **2-hour scan** (recommended): `find_moon_above_horizon` — fast, ~12 evals/day
//! - **10-minute scan** (validation): `find_moon_above_horizon_scan` — finer step
//!
//! Both delegate to `math_core::intervals` for scan + Brent refinement.
//!
//! ## Performance Targets
//!
//! - Single altitude computation: <1ms
//! - 365-day search: <1s

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use qtty::*;
use siderust::bodies::Moon;
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::ObserverSite;
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

    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);

    Period::new(mjd_start, mjd_end)
}

// =============================================================================
// Single Altitude Computation Benchmark
// =============================================================================

fn bench_moon_altitude_computation(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
    let mjd = ModifiedJulianDate::new(51544.5); // J2000

    let mut group = c.benchmark_group("moon_altitude_single");

    // Benchmark single altitude computation via trait
    group.bench_function("compute_altitude", |b| {
        b.iter(|| {
            let _altitude = Moon.altitude_at(black_box(&site), black_box(mjd));
        });
    });

    group.finish();
}

// =============================================================================
// Moon Above Horizon Benchmarks (2-hour scan — recommended)
// =============================================================================

fn bench_moon_above_horizon(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_above_horizon");

    // 1-day horizon
    group.bench_function("find_moon_above_horizon_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let _result = Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    // 7-day horizon
    group.bench_function("find_moon_above_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    // 30-day horizon (full lunar cycle)
    group.bench_function("find_moon_above_horizon_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let _result = Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    // 365-day horizon (full year) - PRIMARY PERFORMANCE TARGET
    group.bench_function("find_moon_above_horizon_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Moon Below Horizon Benchmarks (2-hour scan — recommended)
// =============================================================================

fn bench_moon_below_horizon(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_below_horizon");

    // 1-day horizon
    group.bench_function("find_moon_below_horizon_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let _result = Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
            );
        });
    });

    // 7-day horizon
    group.bench_function("find_moon_below_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
            );
        });
    });

    // 30-day horizon
    group.bench_function("find_moon_below_horizon_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let _result = Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
            );
        });
    });

    // 365-day horizon (full year) - PRIMARY PERFORMANCE TARGET
    group.bench_function("find_moon_below_horizon_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Altitude Range Benchmarks
// =============================================================================

fn bench_moon_altitude_range(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_altitude_range");

    // Finding Moon at low altitude (0-30 degrees) over 7 days
    group.bench_function("find_moon_low_altitude_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = Moon.altitude_periods(&AltitudeQuery {
                observer: black_box(site),
                window: black_box(period),
                min_altitude: black_box(Degrees::new(0.0)),
                max_altitude: black_box(Degrees::new(30.0)),
            });
        });
    });

    // Finding Moon at high altitude (60-90 degrees) over 7 days
    group.bench_function("find_moon_high_altitude_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = Moon.altitude_periods(&AltitudeQuery {
                observer: black_box(site),
                window: black_box(period),
                min_altitude: black_box(Degrees::new(60.0)),
                max_altitude: black_box(Degrees::new(90.0)),
            });
        });
    });

    group.finish();
}

// =============================================================================
// Algorithm Comparison: Cached (default) vs above/below for 365-day horizons
// =============================================================================

fn bench_algorithm_comparison(c: &mut Criterion) {
    let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_algorithm_comparison");
    group.measurement_time(Duration::from_secs(15));

    // 365-day above vs below comparison
    group.bench_function("moon_above_horizon_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = Moon.above_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    group.bench_function("moon_below_horizon_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = Moon.below_threshold(
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    group.finish();
}

criterion_group! {
    name = moon_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(10))
        .sample_size(20);
    targets = bench_moon_altitude_computation,
              bench_moon_above_horizon,
              bench_moon_below_horizon,
              bench_moon_altitude_range,
              bench_algorithm_comparison
}
criterion_main!(moon_benches);
