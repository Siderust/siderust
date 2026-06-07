// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for lunar altitude calculations.
//!
//! Tests the performance of calculating Moon altitude and finding periods when
//! the Moon is above/below the horizon using the ELP2000 lunar theory.
//!
//! ## Algorithms Compared
//!
//! - **2-hour scan** (recommended): `above_threshold` with default options, fast, ~12 evals/day
//! - **10-minute scan** (validation): finer `scan_step_days`, for algorithm comparison
//!
//! Both the default Chebyshev-first engine and the scan+Brent baseline
//! are compared over 30, 184, and 365 day windows.
//!
//! ## Performance Targets
//!
//! - Single altitude computation: <1ms
//! - 365-day search: <1s

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bodies::Moon;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, AltitudeProvider, SearchOpts,
};
use siderust::qtty::*;
use siderust::time::{Interval, ModifiedJulianDate};
use std::hint::black_box;
use std::time::Duration;

fn build_period(days: u32) -> Interval<ModifiedJulianDate> {
    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap())
        + chrono::Duration::days(days as i64);

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    let mjd_start: ModifiedJulianDate = ModifiedJulianDate::from(start_dt);
    let mjd_end: ModifiedJulianDate = ModifiedJulianDate::from(end_dt);

    Interval::new(mjd_start, mjd_end)
}

// =============================================================================
// Single Altitude Computation Benchmark
// =============================================================================

fn bench_moon_altitude_computation(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
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
// Moon Above Horizon Benchmarks (2-hour scan, recommended)
// =============================================================================

fn bench_moon_above_horizon(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();

    let mut group = c.benchmark_group("moon_above_horizon");

    // 1-day horizon
    group.bench_function("above_threshold_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 7-day horizon
    group.bench_function("above_threshold_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 30-day horizon (full lunar cycle)
    group.bench_function("above_threshold_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 365-day horizon (full year) - PRIMARY PERFORMANCE TARGET
    group.bench_function("above_threshold_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Moon Below Horizon Benchmarks (2-hour scan, recommended)
// =============================================================================

fn bench_moon_below_horizon(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();

    let mut group = c.benchmark_group("moon_below_horizon");

    // 1-day horizon
    group.bench_function("below_threshold_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let _result = below_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 7-day horizon
    group.bench_function("below_threshold_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = below_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 30-day horizon
    group.bench_function("below_threshold_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let _result = below_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // 365-day horizon (full year) - PRIMARY PERFORMANCE TARGET
    group.bench_function("below_threshold_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = below_threshold(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(-0.5)),
                black_box(SearchOpts::default()),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Altitude Range Benchmarks
// =============================================================================

fn bench_moon_altitude_range(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();

    let mut group = c.benchmark_group("moon_altitude_range");

    // Finding Moon at low altitude (0-30 degrees) over 7 days
    group.bench_function("altitude_ranges_low_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = altitude_ranges(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(Degrees::new(30.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    // Finding Moon at high altitude (60-90 degrees) over 7 days
    group.bench_function("altitude_ranges_high_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = altitude_ranges(
                black_box(&Moon),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(60.0)),
                black_box(Degrees::new(90.0)),
                black_box(SearchOpts::default()),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Algorithm Comparison: Chebyshev-first vs scan+Brent baseline
// =============================================================================

fn bench_algorithm_comparison(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let default_opts = SearchOpts::default();
    let scan_opts = SearchOpts {
        scan_step_days: Some(Quantity::<Hour>::new(2.0).to::<Day>()),
        ..default_opts
    };

    let mut group = c.benchmark_group("moon_algorithm_comparison");
    group.measurement_time(Duration::from_secs(15));

    for (days_label, days) in [("30day", 30), ("184day", 184), ("365day", 365)] {
        let period = build_period(days);

        for (mode_label, opts) in [("chebyshev_first", default_opts), ("scan_brent", scan_opts)] {
            group.bench_function(format!("above_horizon/{mode_label}/{days_label}"), |b| {
                b.iter(|| {
                    let _result = above_threshold(
                        black_box(&Moon),
                        black_box(&site),
                        black_box(period),
                        black_box(Degrees::new(0.0)),
                        black_box(opts),
                    );
                });
            });

            group.bench_function(format!("below_horizon/{mode_label}/{days_label}"), |b| {
                b.iter(|| {
                    let _result = below_threshold(
                        black_box(&Moon),
                        black_box(&site),
                        black_box(period),
                        black_box(Degrees::new(0.0)),
                        black_box(opts),
                    );
                });
            });
        }
    }

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
