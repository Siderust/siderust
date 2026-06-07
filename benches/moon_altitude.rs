// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks comparing the lunar Chebyshev context engine against scan+Brent baseline.
//!
//! Requires the `bench-internals` feature:
//! `cargo bench --features bench-internals --bench moon_altitude`

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bench_internals;
use siderust::bodies::Moon;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{above_threshold, below_threshold, SearchOpts};
use siderust::qtty::Degrees;
use siderust::time::{Interval, ModifiedJulianDate};
use std::hint::black_box;
use std::time::Duration;

fn build_period(days: u32) -> Interval<ModifiedJulianDate> {
    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = start_naive + chrono::Duration::days(days as i64);
    Interval::new(
        ModifiedJulianDate::from(Utc.from_utc_datetime(&start_naive)),
        ModifiedJulianDate::from(Utc.from_utc_datetime(&end_naive)),
    )
}

fn bench_moon_engines(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        let mut cheb = c.benchmark_group(format!("moon/chebyshev_context/{label}"));
        cheb.bench_function("above_threshold/horizon", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ));
            });
        });
        cheb.bench_function("below_threshold/horizon", |b| {
            b.iter(|| {
                black_box(below_threshold(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ));
            });
        });
        cheb.finish();

        let mut scan = c.benchmark_group(format!("moon/scan_brent_baseline/{label}"));
        scan.bench_function("above_threshold/horizon", |b| {
            b.iter(|| {
                black_box(bench_internals::lunar_above_threshold_scan_baseline(
                    site,
                    black_box(period),
                    Degrees::new(0.0),
                    opts,
                ));
            });
        });
        scan.finish();
    }
}

criterion_group! {
    name = moon_benches;
    config = Criterion::default().measurement_time(Duration::from_secs(5)).sample_size(20);
    targets = bench_moon_engines
}
criterion_main!(moon_benches);
