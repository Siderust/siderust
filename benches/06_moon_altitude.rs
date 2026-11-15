// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lunar altitude benchmarks.
//!
//! Covers cached evaluations, short event windows and long search-engine comparisons.
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
//!     window using the same interpolated altitude context
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
use siderust::astro::nutation::NutationModel;
use siderust::bench_internals::{self, MoonAltitudeContext, MoonPositionCache};
use siderust::bodies::Moon;
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::{
    centers::Geodetic,
    frames::{EclipticMeanJ2000, ECEF},
    transform::{context::DefaultNutationModel, frame_rotation_with, AstroContext},
};
use siderust::event::altitude::{
    above_threshold, altitude_ranges, below_threshold, crossings, SearchOpts,
};
use siderust::qtty::{Degrees, DEG, M};
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

// January 2024 in TT, within compact ephemeris coverage.
fn audit_period(days: f64) -> Interval<ModifiedJulianDate> {
    Interval::new(
        ModifiedJulianDate::new(60_310.0),
        ModifiedJulianDate::new(60_310.0 + days),
    )
}

fn bench_moon_evaluation(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let period = audit_period(30.0);
    let cache = MoonAltitudeContext::new(period.start, period.end, site);
    let positions = MoonPositionCache::new(period.start, period.end);
    let context: AstroContext = AstroContext::default();
    let times: Vec<_> = (0..1024)
        .map(|i| ModifiedJulianDate::new(period.start.value() + i as f64 * 30.0 / 1024.0))
        .collect();
    let epochs: Vec<_> = times.iter().map(|t| t.to::<siderust::time::JD>()).collect();

    let mut group = c.benchmark_group("evaluation/moon");
    group.bench_function("cached_altitude", |b| {
        let mut index = 0;
        b.iter(|| {
            let mjd = times[index % times.len()];
            index += 1;
            black_box(cache.altitude_rad(black_box(mjd)))
        })
    });
    group.bench_function("cached_position", |b| {
        let mut index = 0;
        b.iter(|| {
            let mjd = times[index % times.len()];
            index += 1;
            black_box(positions.get_position_km(black_box(mjd)))
        })
    });
    group.bench_function("earth_orientation", |b| {
        let mut index = 0;
        b.iter(|| {
            let jd = epochs[index % epochs.len()];
            index += 1;
            black_box(frame_rotation_with::<EclipticMeanJ2000, ECEF, _>(
                black_box(jd),
                &context,
            ))
        })
    });
    group.bench_function("eop", |b| {
        let mut index = 0;
        b.iter(|| {
            let jd = epochs[index % epochs.len()];
            index += 1;
            black_box(context.eop_at_tt(black_box(jd)))
        })
    });
    group.bench_function("nutation", |b| {
        let mut index = 0;
        b.iter(|| {
            let jd = epochs[index % epochs.len()];
            index += 1;
            black_box(DefaultNutationModel::nutation(black_box(jd)))
        })
    });
    for days in [1.0, 7.0, 30.0] {
        let window = audit_period(days);
        group.bench_function(format!("build_context/{days}d"), |b| {
            b.iter(|| {
                black_box(MoonAltitudeContext::new(
                    black_box(window.start),
                    black_box(window.end),
                    site,
                ))
            })
        });
    }
    group.finish();
}

fn bench_moon_event_windows(c: &mut Criterion) {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let polar = Geodetic::<ECEF>::new(18.95 * DEG, 69.65 * DEG, 0.0 * M);
    let opts = SearchOpts::default();
    for days in [1.0, 7.0, 30.0] {
        let window = audit_period(days);
        let mut group = c.benchmark_group(format!("event_windows/moon/{days}d"));
        group.bench_function("rise_set", |b| {
            b.iter(|| black_box(crossings(&Moon, &site, black_box(window), 0.0 * DEG, opts)))
        });
        group.bench_function("above_threshold", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Moon,
                    &site,
                    black_box(window),
                    0.0 * DEG,
                    opts,
                ))
            })
        });
        group.bench_function("altitude_band", |b| {
            b.iter(|| {
                black_box(altitude_ranges(
                    &Moon,
                    &site,
                    black_box(window),
                    0.0 * DEG,
                    30.0 * DEG,
                    opts,
                ))
            })
        });
        group.bench_function("high_latitude_above_threshold", |b| {
            b.iter(|| {
                black_box(above_threshold(
                    &Moon,
                    &polar,
                    black_box(window),
                    0.0 * DEG,
                    opts,
                ))
            })
        });
        group.finish();
    }

    // The Moon briefly dips below the horizon near lower culmination.
    let shallow = Interval::new(
        ModifiedJulianDate::new(60_335.4),
        ModifiedJulianDate::new(60_335.6),
    );
    c.bench_function("event_windows/moon/shallow_horizon", |b| {
        b.iter(|| {
            black_box(crossings(
                &Moon,
                &polar,
                black_box(shallow),
                0.0 * DEG,
                opts,
            ))
        })
    });
}

criterion_group! {
    name = moon_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);
    targets = bench_moon_evaluation, bench_moon_event_windows, bench_public_api_moon, bench_engines_moon
}
criterion_main!(moon_benches);
