// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Full ephemeris queries and optional compact position comparisons.
//!
//! Run with:
//! ```bash
//! cargo bench --bench ephemeris_comparison
//! ```

use criterion::{criterion_group, criterion_main, Criterion};
#[cfg(feature = "embedded-compact")]
use siderust::ephemeris::CompactEphemeris;
use siderust::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::qtty::Days;
use std::hint::black_box;
use std::time::Duration;

// =============================================================================
// Sun Barycentric
// =============================================================================

fn bench_sun_barycentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/sun_barycentric");

    let jd = siderust::time::J2000 + Days::new(25.0 * 365.25);
    group.bench_function("vsop87", |b| {
        b.iter(|| Vsop87Ephemeris::sun_barycentric(black_box(jd)));
    });
    #[cfg(feature = "embedded-compact")]
    group.bench_function("compact", |b| {
        CompactEphemeris::try_sun_barycentric(jd)
            .expect("benchmark epoch must be in compact coverage");
        b.iter(|| CompactEphemeris::sun_barycentric(black_box(jd)));
    });

    group.finish();
}

// =============================================================================
// Earth Heliocentric
// =============================================================================

fn bench_earth_heliocentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/earth_heliocentric");

    let jd = siderust::time::J2000 + Days::new(25.0 * 365.25);
    group.bench_function("vsop87", |b| {
        b.iter(|| Vsop87Ephemeris::earth_heliocentric(black_box(jd)));
    });
    #[cfg(feature = "embedded-compact")]
    group.bench_function("compact", |b| {
        CompactEphemeris::try_earth_heliocentric(jd)
            .expect("benchmark epoch must be in compact coverage");
        b.iter(|| CompactEphemeris::earth_heliocentric(black_box(jd)));
    });

    group.finish();
}

// =============================================================================
// Earth Barycentric Velocity
// =============================================================================

fn bench_earth_velocity(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/earth_velocity");

    let jd = siderust::time::J2000 + Days::new(25.0 * 365.25);
    group.bench_function("vsop87", |b| {
        b.iter(|| Vsop87Ephemeris::earth_barycentric_velocity(black_box(jd)));
    });

    group.finish();
}

// =============================================================================
// Moon Geocentric
// =============================================================================

fn bench_moon_geocentric(c: &mut Criterion) {
    let mut group = c.benchmark_group("ephemeris/moon_geocentric");

    let jd = siderust::time::J2000 + Days::new(25.0 * 365.25);
    group.bench_function("vsop87", |b| {
        b.iter(|| Vsop87Ephemeris::moon_geocentric(black_box(jd)));
    });
    #[cfg(feature = "embedded-compact")]
    group.bench_function("compact", |b| {
        CompactEphemeris::try_moon_geocentric(jd)
            .expect("benchmark epoch must be in compact coverage");
        b.iter(|| CompactEphemeris::moon_geocentric(black_box(jd)));
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
