// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks comparing the solar daily predictor against internal baselines.

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bench_internals;
use siderust::bodies::Sun;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{altitude_ranges, below_threshold, crossings, SearchOpts};
use siderust::event::solar::twilight;
use siderust::qtty::Degrees;
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

    Interval::new(
        ModifiedJulianDate::from(start_dt),
        ModifiedJulianDate::from(end_dt),
    )
}

fn bench_solar_engines(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        let mut daily = c.benchmark_group(format!("solar/daily_predictor/{label}"));
        daily.bench_function("below_threshold/-18", |b| {
            b.iter(|| {
                black_box(below_threshold(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(twilight::ASTRONOMICAL),
                    opts,
                ));
            });
        });
        daily.bench_function("altitude_ranges/twilight", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(-18.0)),
                    black_box(Degrees::new(-12.0)),
                    opts,
                ));
            });
        });
        daily.bench_function("crossings/horizon", |b| {
            b.iter(|| {
                black_box(crossings(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ));
            });
        });
        daily.finish();

        let mut cheb = c.benchmark_group(format!("solar/chebyshev_fallback_baseline/{label}"));
        cheb.bench_function("below_threshold/-18", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_below_threshold_chebyshev_baseline(
                    site,
                    black_box(period),
                    twilight::ASTRONOMICAL,
                    opts,
                ));
            });
        });
        cheb.finish();

        let mut scan = c.benchmark_group(format!("solar/scan_brent_baseline/{label}"));
        scan.bench_function("below_threshold/-18", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_below_threshold_scan_baseline(
                    site,
                    black_box(period),
                    twilight::ASTRONOMICAL,
                    opts,
                ));
            });
        });
        scan.finish();
    }
}

criterion_group! {
    name = solar_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_solar_engines
}
criterion_main!(solar_benches);
