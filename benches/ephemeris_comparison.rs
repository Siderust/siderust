// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for the `Vsop87Ephemeris` backend.
//!
//! Run with:
//! ```bash
//! cargo bench --bench ephemeris_comparison
//! ```

use criterion::{criterion_group, criterion_main, Criterion};
use siderust::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::qtty::Days;
use std::hint::black_box;
use std::time::Duration;

// =============================================================================
// Sun Barycentric
// =============================================================================

fn bench_sun_barycentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/sun_barycentric");

    group.bench_function("vsop87", |b| {
        let mut jd = siderust::time::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::sun_barycentric(black_box(jd));
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
        let mut jd = siderust::time::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::earth_heliocentric(black_box(jd));
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
        let mut jd = siderust::time::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::earth_barycentric_velocity(black_box(jd));
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
        let mut jd = siderust::time::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = Vsop87Ephemeris::moon_geocentric(black_box(jd));
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
