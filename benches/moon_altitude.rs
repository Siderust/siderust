// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lunar altitude benchmarks.
//!
//! Organized into two halves:
//!
//! - **`public_api/moon/{30d,184d,365d}`**: what callers invoke via the stable
//!   Option A API. The default engine (`MoonAltitudeContext` + Chebyshev-first
//!   crossing discovery + local scan+Brent fallback per unsafe segment) runs here.
//!
//! - **`engines/moon/{30d,184d,365d}`**: apples-to-apples internal engine
//!   comparison using `bench_internals` baselines:
//!   - `chebyshev_context` — default engine (same as `public_api`, labelled for
//!     direct comparison)
//!   - `scan_brent_baseline` — uniform scan + Brent refinement over the full
//!     window (not exposed publicly)
//!
//! Note: topocentric parallax (~1° at the horizon) means each precise lunar
//! altitude evaluation is heavier than the solar equivalent.
//!
//! Requires the `bench-internals` feature:
//!
//! ```bash
//! cargo bench --features bench-internals --bench moon_altitude
//! cargo bench --no-run --features bench-internals
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bench_internals;
use siderust::bodies::Moon;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::event::altitude::{above_threshold, altitude_ranges, below_threshold, SearchOpts};
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

fn bench_public_api_moon(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        let mut g = c.benchmark_group(format!("public_api/moon/{label}"));

        g.bench_function("above_threshold/horizon", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ))
            });
        });

        g.bench_function("below_threshold/horizon", |b| {
            b.iter(|| {
                black_box(below_threshold(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ))
            });
        });

        g.bench_function("altitude_ranges/observation_window", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    black_box(Degrees::new(30.0)),
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

fn bench_engines_moon(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        // --- above_threshold: Chebyshev context vs scan+Brent ------------------
        let mut g = c.benchmark_group(format!("engines/moon/{label}/above_threshold"));

        // Default engine (matches public_api group above)
        g.bench_function("chebyshev_context", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    opts,
                ))
            });
        });

        g.bench_function("scan_brent_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::lunar_above_threshold_scan_baseline(
                    site,
                    black_box(period),
                    Degrees::new(0.0),
                    opts,
                ))
            });
        });

        g.finish();

        // --- altitude_ranges: Chebyshev context vs scan+Brent ------------------
        let mut g = c.benchmark_group(format!("engines/moon/{label}/altitude_ranges"));

        g.bench_function("chebyshev_context", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Moon,
                    &site,
                    black_box(period),
                    black_box(Degrees::new(0.0)),
                    black_box(Degrees::new(30.0)),
                    opts,
                ))
            });
        });

        g.bench_function("scan_brent_baseline", |b| {
            b.iter(|| {
                black_box(bench_internals::lunar_altitude_ranges_scan_baseline(
                    site,
                    black_box(period),
                    Degrees::new(0.0),
                    Degrees::new(30.0),
                    opts,
                ))
            });
        });

        g.finish();
    }
}

criterion_group! {
    name = moon_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_public_api_moon, bench_engines_moon
}
criterion_main!(moon_benches);
