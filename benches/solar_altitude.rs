// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for solar altitude period calculations.
//!
//! Tests the performance of finding astronomical night periods using different
//! algorithms and time horizons.

use criterion::{criterion_group, criterion_main, Criterion};
use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use siderust::calculus::solar::altitude_periods::{
    find_day_periods, find_night_periods, find_night_periods_scan,
};
use siderust::calculus::solar::twilight;
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

fn bench_find_night_periods(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("solar_altitude_periods");

    // Benchmark for 1-day horizon
    group.bench_function("find_night_periods_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let _result = find_night_periods(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    // Benchmark for 7-day horizon
    group.bench_function("find_night_periods_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = find_night_periods(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    // Benchmark for 30-day horizon
    group.bench_function("find_night_periods_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let _result = find_night_periods(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    // Benchmark for 365-day horizon (full year)
    group.bench_function("find_night_periods_365day", |b| {
        let period = black_box(build_period(365));
        b.iter(|| {
            let _result = find_night_periods(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    // Compare with scan-based algorithm for 7 days
    group.bench_function("find_night_periods_scan_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = find_night_periods_scan(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    // Compare with culmination-based algorithm for 7 days
    group.bench_function("find_day_periods_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let _result = find_day_periods(
                black_box(site),
                black_box(period),
                black_box(twilight::ASTRONOMICAL),
            );
        });
    });

    group.finish();
}

criterion_group! {
    name = solar_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_find_night_periods
}
criterion_main!(solar_benches);
