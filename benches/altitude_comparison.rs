// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cross-body altitude comparison benchmarks.
//!
//! Measures single-point altitude evaluation and period-finding performance for
//! Sun, Moon, and Star (Sirius) side-by-side. Criterion's HTML reports show the
//! relative cost of each body's altitude engine.
//!
//! Group names follow the project taxonomy:
//!
//! - `altitude/single_eval` — single call to `altitude_at`
//! - `public_api/comparison/above_horizon_7d` — 7-day above-horizon search
//! - `public_api/comparison/below_threshold_30d` — 30-day below-threshold
//! - `public_api/comparison/above_horizon_365d` — full-year above-horizon
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench altitude_comparison
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bodies::{Moon, Sun};
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::spherical::direction;
use siderust::event::altitude::{above_threshold, below_threshold, AltitudeProvider, SearchOpts};
use siderust::qtty::*;
use siderust::time::{Interval, ModifiedJulianDate};
use std::hint::black_box;
use std::time::Duration;

fn build_period(days: u32) -> Interval<ModifiedJulianDate> {
    let start = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end = start + chrono::Duration::days(days as i64);
    Interval::new(
        ModifiedJulianDate::from(Utc.from_utc_datetime(&start)),
        ModifiedJulianDate::from(Utc.from_utc_datetime(&end)),
    )
}

fn sirius_icrs() -> direction::ICRS {
    direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716))
}

// ---------------------------------------------------------------------------
// Single-point altitude evaluation
// ---------------------------------------------------------------------------

fn bench_single_altitude(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let mjd = ModifiedJulianDate::new(51544.5); // J2000
    let sirius = sirius_icrs();

    let mut g = c.benchmark_group("altitude/single_eval");

    g.bench_function("sun", |b| {
        b.iter(|| Sun.altitude_at(black_box(&site), black_box(mjd)));
    });
    g.bench_function("moon", |b| {
        b.iter(|| Moon.altitude_at(black_box(&site), black_box(mjd)));
    });
    g.bench_function("star_sirius", |b| {
        b.iter(|| sirius.altitude_at(black_box(&site), black_box(mjd)));
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 7-day above-horizon search
// ---------------------------------------------------------------------------

fn bench_above_horizon_7d(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let period = build_period(7);
    let sirius = sirius_icrs();
    let opts = SearchOpts::default();

    let mut g = c.benchmark_group("public_api/comparison/above_horizon_7d");

    g.bench_function("sun", |b| {
        b.iter(|| {
            above_threshold(
                &Sun,
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });
    g.bench_function("moon", |b| {
        b.iter(|| {
            above_threshold(
                &Moon,
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });
    g.bench_function("star_sirius", |b| {
        b.iter(|| {
            above_threshold(
                black_box(&sirius),
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 30-day below-threshold search
// ---------------------------------------------------------------------------

fn bench_below_threshold_30d(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let period = build_period(30);
    let sirius = sirius_icrs();
    let opts = SearchOpts::default();

    let mut g = c.benchmark_group("public_api/comparison/below_threshold_30d");

    g.bench_function("sun_astro_night", |b| {
        b.iter(|| {
            below_threshold(
                &Sun,
                &site,
                black_box(period),
                black_box(Degrees::new(-18.0)),
                opts,
            )
        });
    });
    g.bench_function("moon_below_horizon", |b| {
        b.iter(|| {
            below_threshold(
                &Moon,
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });
    g.bench_function("star_sirius_below_horizon", |b| {
        b.iter(|| {
            below_threshold(
                black_box(&sirius),
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 365-day above-horizon search (primary performance target)
// ---------------------------------------------------------------------------

fn bench_above_horizon_365d(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let period = build_period(365);
    let sirius = sirius_icrs();
    let opts = SearchOpts::default();

    let mut g = c.benchmark_group("public_api/comparison/above_horizon_365d");

    g.bench_function("sun", |b| {
        b.iter(|| {
            above_threshold(
                &Sun,
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });
    g.bench_function("moon", |b| {
        b.iter(|| {
            above_threshold(
                &Moon,
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });
    g.bench_function("star_sirius", |b| {
        b.iter(|| {
            above_threshold(
                black_box(&sirius),
                &site,
                black_box(period),
                black_box(Degrees::new(0.0)),
                opts,
            )
        });
    });

    g.finish();
}

criterion_group! {
    name = altitude_comparison_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(10))
        .sample_size(20)
        .without_plots();
    targets = bench_single_altitude,
              bench_above_horizon_7d,
              bench_below_threshold_30d,
              bench_above_horizon_365d
}
criterion_main!(altitude_comparison_benches);
