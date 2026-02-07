// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for fixed-star altitude calculations.
//!
//! Uses the unified altitude API with a static J2000 RA/Dec target
//! (Sirius) to measure single‑point altitude and range/threshold searches.

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use qtty::*;
use siderust::calculus::altitude::{
    above_threshold, altitude_at, altitude_ranges, crossings, AltitudeTarget, SearchOpts,
};
use siderust::calculus::stellar;
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};
use std::hint::black_box;
use std::time::Duration;

fn build_period(days: u32) -> Period<ModifiedJulianDate> {
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

fn sirius_target() -> AltitudeTarget {
    AltitudeTarget::FixedEquatorial {
        ra: Degrees::new(101.287),
        dec: Degrees::new(-16.716),
    }
}

// =============================================================================
// Single Altitude Computation Benchmark
// =============================================================================

fn bench_star_altitude_computation(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    let target = sirius_target();
    let jd = siderust::astro::JulianDate::J2000;

    let mut group = c.benchmark_group("star_altitude_single");

    group.bench_function("compute_altitude", |b| {
        b.iter(|| {
            let _altitude = altitude_at(black_box(&target), black_box(&site), black_box(jd));
        });
    });

    group.finish();
}

// =============================================================================
// Threshold + Range Benchmarks
// =============================================================================

fn bench_star_thresholds(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    let target = sirius_target();
    let opts = SearchOpts::default();

    let mut group = c.benchmark_group("star_altitude_thresholds");

    group.bench_function("above_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&target),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(opts),
            );
        });
    });

    group.bench_function("above_30deg_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = above_threshold(
                black_box(&target),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(30.0)),
                black_box(opts),
            );
        });
    });

    group.bench_function("altitude_range_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = altitude_ranges(
                black_box(&target),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(10.0)),
                black_box(Degrees::new(50.0)),
                black_box(opts),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Crossing Benchmarks
// =============================================================================

fn bench_star_crossings(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    let target = sirius_target();
    let opts = SearchOpts::default();

    let mut group = c.benchmark_group("star_altitude_crossings");

    group.bench_function("crossings_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = crossings(
                black_box(&target),
                black_box(&site),
                black_box(period),
                black_box(Degrees::new(0.0)),
                black_box(opts),
            );
        });
    });

    group.finish();
}

// =============================================================================
// Analytical vs Scan Comparison
// =============================================================================

fn bench_stellar_analytical_vs_scan(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    let ra = Degrees::new(101.287);
    let dec = Degrees::new(-16.716);

    let mut group = c.benchmark_group("stellar_analytical_vs_scan");

    group.bench_function("analytical_above_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = stellar::find_star_above_periods(
                black_box(ra),
                black_box(dec),
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    group.bench_function("scan_above_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = stellar::find_star_above_periods_scan(
                black_box(ra),
                black_box(dec),
                black_box(site),
                black_box(period),
                black_box(Degrees::new(0.0)),
            );
        });
    });

    group.bench_function("analytical_range_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = stellar::find_star_range_periods(
                black_box(ra),
                black_box(dec),
                black_box(site),
                black_box(period),
                black_box((Degrees::new(10.0), Degrees::new(50.0))),
            );
        });
    });

    group.bench_function("scan_range_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = stellar::find_star_range_periods_scan(
                black_box(ra),
                black_box(dec),
                black_box(site),
                black_box(period),
                black_box((Degrees::new(10.0), Degrees::new(50.0))),
            );
        });
    });

    group.finish();
}

criterion_group! {
    name = star_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_star_altitude_computation, bench_star_thresholds, bench_star_crossings, bench_stellar_analytical_vs_scan
}
criterion_main!(star_benches);
