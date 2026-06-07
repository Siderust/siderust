// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for solar altitude period calculations.
//!
//! Compares the default Chebyshev-first engine against the legacy scan+Brent
//! baseline over 30, 184, and 365 day windows.

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use siderust::bodies::Sun;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{
    altitude_ranges, below_threshold, AltitudePeriodsProvider, AltitudeQuery, SearchOpts,
};
use siderust::event::solar::twilight;
use siderust::qtty::{Day, Degrees, Hour, Quantity};
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

fn bench_find_night_periods(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let default_opts = SearchOpts::default();
    let scan_opts = SearchOpts {
        scan_step_days: Some(Quantity::<Hour>::new(2.0).to::<Day>()),
        ..default_opts
    };

    let mut group = c.benchmark_group("solar_altitude_periods");

    for (label, days) in [("1month", 30), ("6months", 184), ("1year", 365)] {
        let period = build_period(days);

        group.bench_function(
            BenchmarkId::new("below_threshold/chebyshev_first", label),
            |b| {
                b.iter(|| {
                    let _result = below_threshold(
                        black_box(&Sun),
                        black_box(&site),
                        black_box(period),
                        black_box(twilight::ASTRONOMICAL),
                        black_box(default_opts),
                    );
                });
            },
        );

        group.bench_function(BenchmarkId::new("below_threshold/scan_brent", label), |b| {
            b.iter(|| {
                let _result = below_threshold(
                    black_box(&Sun),
                    black_box(&site),
                    black_box(period),
                    black_box(twilight::ASTRONOMICAL),
                    black_box(scan_opts),
                );
            });
        });

        group.bench_function(
            BenchmarkId::new("altitude_periods/chebyshev_first", label),
            |b| {
                b.iter(|| {
                    let _result = altitude_ranges(
                        black_box(&Sun),
                        black_box(&site),
                        black_box(period),
                        black_box(Degrees::new(-90.0)),
                        black_box(twilight::ASTRONOMICAL),
                        black_box(default_opts),
                    );
                });
            },
        );

        group.bench_function(
            BenchmarkId::new("altitude_periods/scan_brent", label),
            |b| {
                b.iter(|| {
                    let _result = altitude_ranges(
                        black_box(&Sun),
                        black_box(&site),
                        black_box(period),
                        black_box(Degrees::new(-90.0)),
                        black_box(twilight::ASTRONOMICAL),
                        black_box(scan_opts),
                    );
                });
            },
        );

        group.bench_function(BenchmarkId::new("provider_default", label), |b| {
            b.iter(|| {
                let query = AltitudeQuery {
                    observer: black_box(site),
                    window: black_box(period),
                    min_altitude: black_box(Degrees::new(-90.0)),
                    max_altitude: black_box(twilight::ASTRONOMICAL),
                    correction_policy: siderust::astro::apparent::CorrectionPolicy::APPARENT,
                };
                let _result = Sun.altitude_periods(black_box(&query));
            });
        });
    }

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
