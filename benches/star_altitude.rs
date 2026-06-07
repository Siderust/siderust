// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Fixed-star altitude benchmarks.
//!
//! Organized as:
//!
//! - **`public_api/star/{30d,184d,365d}`**: stable Option A API calls for a
//!   fixed ICRS target (Sirius, α CMa). The star engine uses an analytic
//!   sinusoidal model exploiting Earth's diurnal rotation; no internal baseline
//!   variants exist.
//!
//! - **`altitude/single_eval`**: single-point altitude cost for the star engine,
//!   useful as a reference for the period-finding benchmarks above.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench star_altitude
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::spherical::direction;
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, crossings, AltitudeProvider, SearchOpts,
};
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

/// Sirius (α CMa): RA 101.287°, Dec −16.716°.
fn sirius() -> direction::ICRS {
    direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716))
}

// ---------------------------------------------------------------------------
// Single-point altitude evaluation (reference baseline)
// ---------------------------------------------------------------------------

fn bench_star_single_eval(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let target = sirius();
    let mjd = ModifiedJulianDate::new(51544.5); // J2000

    let mut g = c.benchmark_group("altitude/single_eval/star");
    g.bench_function("sirius", |b| {
        b.iter(|| target.altitude_at(black_box(&site), black_box(mjd)));
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// A. Public API benchmarks
// ---------------------------------------------------------------------------

fn bench_public_api_star(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let target = sirius();
    let opts = SearchOpts::default();

    for (label, days) in [("30d", 30u32), ("184d", 184), ("365d", 365)] {
        let period = build_period(days);

        let mut g = c.benchmark_group(format!("public_api/star/{label}"));

        g.bench_function("above_threshold/horizon", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    black_box(&target),
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
                    black_box(&target),
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
                    black_box(&target),
                    &site,
                    black_box(period),
                    black_box(Degrees::new(10.0)),
                    black_box(Degrees::new(50.0)),
                    opts,
                ))
            });
        });

        g.bench_function("crossings/horizon", |b| {
            b.iter(|| {
                black_box(crossings(
                    black_box(&target),
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

criterion_group! {
    name = star_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_star_single_eval, bench_public_api_star
}
criterion_main!(star_benches);
