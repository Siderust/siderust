// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Comparative benchmarks across ephemeris backends.
//!
//! Measures every `Ephemeris` trait method for VSOP87, DE440 (if enabled), and
//! DE441 (if enabled) side-by-side in the same Criterion group so that
//! Criterion's HTML reports show relative performance.
//!
//! Run with:
//! ```bash
//! cargo bench --bench ephemeris_comparison
//! cargo bench --bench ephemeris_comparison --features de440
//! cargo bench --bench ephemeris_comparison --features de441
//! cargo bench --bench ephemeris_comparison --features de440,de441
//! ```

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::Days;
use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::time::JulianDate;
use std::hint::black_box;
use std::time::Duration;

#[cfg(feature = "de440")]
use siderust::calculus::ephemeris::De440Ephemeris;
#[cfg(feature = "de441")]
use siderust::calculus::ephemeris::De441Ephemeris;

// =============================================================================
// Sun Barycentric
// =============================================================================

fn bench_sun_barycentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/sun_barycentric");

    group.bench_function("vsop87", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::sun_barycentric(black_box(jd));
        });
    });

    #[cfg(feature = "de440")]
    group.bench_function("de440", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De440Ephemeris::sun_barycentric(black_box(jd));
        });
    });

    #[cfg(feature = "de441")]
    group.bench_function("de441", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::sun_barycentric(black_box(jd));
        });
    });

    group.finish();
}

// =============================================================================
// Earth Heliocentric
// =============================================================================

fn bench_earth_heliocentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/earth_heliocentric");

    group.bench_function("vsop87", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::earth_heliocentric(black_box(jd));
        });
    });

    #[cfg(feature = "de440")]
    group.bench_function("de440", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De440Ephemeris::earth_heliocentric(black_box(jd));
        });
    });

    #[cfg(feature = "de441")]
    group.bench_function("de441", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::earth_heliocentric(black_box(jd));
        });
    });

    group.finish();
}

// =============================================================================
// Earth Barycentric Velocity
// =============================================================================

fn bench_earth_velocity(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/earth_velocity");

    group.bench_function("vsop87", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::earth_barycentric_velocity(black_box(jd));
        });
    });

    #[cfg(feature = "de440")]
    group.bench_function("de440", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De440Ephemeris::earth_barycentric_velocity(black_box(jd));
        });
    });

    #[cfg(feature = "de441")]
    group.bench_function("de441", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::earth_barycentric_velocity(black_box(jd));
        });
    });

    group.finish();
}

// =============================================================================
// Moon Geocentric
// =============================================================================

fn bench_moon_geocentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/moon_geocentric");

    group.bench_function("vsop87", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::moon_geocentric(black_box(jd));
        });
    });

    #[cfg(feature = "de440")]
    group.bench_function("de440", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De440Ephemeris::moon_geocentric(black_box(jd));
        });
    });

    #[cfg(feature = "de441")]
    group.bench_function("de441", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::moon_geocentric(black_box(jd));
        });
    });

    group.finish();
}

criterion_group! {
    name = ephemeris_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(5))
        .sample_size(100)
        .without_plots();
    targets = bench_sun_barycentric,
              bench_earth_heliocentric,
              bench_earth_velocity,
              bench_moon_geocentric
}
criterion_main!(ephemeris_benches);
