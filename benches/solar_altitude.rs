// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Solar altitude benchmarks.
//!
//! Organized into two halves:
//!
//! - **`public_api/solar/{30d,184d,365d}`**: what callers invoke via the stable
//!   Option A API. All four public altitude functions are covered. The default
//!   engine (solar daily predictor + precise Brent validation + batch range path)
//!   is what runs here.
//!
//! - **`engines/solar/{30d,184d,365d}`**: apples-to-apples internal engine
//!   comparison for the same operations, using `bench_internals` baselines:
//!   - `daily_predictor` — default engine (same as `public_api`, labelled for
//!     direct comparison)
//!   - `chebyshev_baseline` — generic Chebyshev-first engine, daily predictor
//!     disabled (not exposed publicly)
//!   - `scan_brent_baseline` — uniform scan + Brent refinement (not exposed
//!     publicly)
//!
//! Requires the `bench-internals` feature:
//!
//! ```bash
//! cargo bench --features bench-internals --bench solar_altitude
//! cargo bench --no-run --features bench-internals
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bench_internals;
use siderust::bodies::Sun;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, crossings, SearchOpts,
};
use siderust::event::solar::twilight;
use siderust::qtty::Degrees;
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

// ---------------------------------------------------------------------------
// A. Public API benchmarks — stable Option A calls, default engine
// ---------------------------------------------------------------------------

fn bench_public_api_solar(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        let mut g = c.benchmark_group(format!("public_api/solar/{label}"));

        g.bench_function("below_threshold/astro_night", |b| {
            b.iter(|| {
                black_box(below_threshold(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(twilight::ASTRONOMICAL),
                    opts,
                ))
            });
        });

        g.bench_function("above_threshold/horizon", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ))
            });
        });

        // altitude_ranges uses the batch daily-predictor path internally
        g.bench_function("altitude_ranges/astro_to_nautical", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(-18.0)),
                    black_box(Degrees::new(-12.0)),
                    opts,
                ))
            });
        });

        g.bench_function("crossings/horizon", |b| {
            b.iter(|| {
                black_box(crossings(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ))
            });
        });

        g.finish();
    }
}

// ---------------------------------------------------------------------------
// B. Engine comparison — internal baselines via bench_internals
// ---------------------------------------------------------------------------

fn bench_engines_solar(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        // --- below_threshold: daily predictor vs Chebyshev vs scan+Brent --------
        let mut g = c.benchmark_group(format!("engines/solar/{label}/below_threshold"));

        // Default engine (matches public_api group above)
        g.bench_function("daily_predictor", |b| {
            b.iter(|| {
                black_box(below_threshold(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(twilight::ASTRONOMICAL),
                    opts,
                ))
            });
        });

        g.bench_function("chebyshev_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_below_threshold_chebyshev_baseline(
                    site,
                    black_box(period),
                    twilight::ASTRONOMICAL,
                    opts,
                ))
            });
        });

        g.bench_function("scan_brent_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_below_threshold_scan_baseline(
                    site,
                    black_box(period),
                    twilight::ASTRONOMICAL,
                    opts,
                ))
            });
        });

        g.finish();

        // --- altitude_ranges: batch daily path vs Chebyshev vs scan+Brent -------
        let mut g = c.benchmark_group(format!("engines/solar/{label}/altitude_ranges"));

        // Default engine — uses the batch solar daily path (two thresholds, one pass)
        g.bench_function("daily_batch", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Sun,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(-18.0)),
                    black_box(Degrees::new(-12.0)),
                    opts,
                ))
            });
        });

        g.bench_function("chebyshev_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_altitude_ranges_chebyshev_baseline(
                    site,
                    black_box(period),
                    Degrees::new(-18.0),
                    Degrees::new(-12.0),
                    opts,
                ))
            });
        });

        g.bench_function("scan_brent_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_altitude_ranges_scan_baseline(
                    site,
                    black_box(period),
                    Degrees::new(-18.0),
                    Degrees::new(-12.0),
                    opts,
                ))
            });
        });

        g.finish();
    }
}

// ---------------------------------------------------------------------------
// C. Twilight profile — batch vs independent: demonstrates one-pass advantage
// ---------------------------------------------------------------------------

fn bench_twilight_profile(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    // All five standard twilight thresholds processed in one daily pass.
    let thresholds: &[Degrees] = &[
        Degrees::new(0.0),
        Degrees::new(-0.833),
        Degrees::new(-6.0),
        Degrees::new(-12.0),
        Degrees::new(-18.0),
    ];

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);
        let mut g = c.benchmark_group(format!("engines/solar/{label}/twilight_profile"));

        // Single-pass batch: computes all five thresholds in one daily sweep.
        g.bench_function("daily_batch_all_thresholds", |b| {
            b.iter(|| {
                black_box(bench_internals::solar_twilight_profile(
                    site,
                    black_box(period),
                    thresholds,
                    opts,
                ))
            });
        });

        // Independent baseline: five separate `below_threshold` calls.
        g.bench_function("independent_below_threshold_x5", |b| {
            b.iter(|| {
                black_box(
                    thresholds
                        .iter()
                        .map(|&thr| below_threshold(&Sun, &site, black_box(period), thr, opts))
                        .collect::<Vec<_>>(),
                )
            });
        });

        g.finish();
    }
}

criterion_group! {
    name = solar_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_public_api_solar, bench_engines_solar, bench_twilight_profile
}
criterion_main!(solar_benches);
