// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Micro-benchmarks for ELP2000 evaluation.
//!
//! Run with: `cargo bench --bench elp2000`
//!
//! Baseline (before SIMD optimization):
//! - get_geo_position_j2000: ~937 µs
//! - get_geo_position_offset: ~1.66 ms

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::{Days, Kilometer};
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::Moon;
use std::hint::black_box;

fn bench_moon_geo_position(c: &mut Criterion) {
    let mut group = c.benchmark_group("elp2000_moon_geo");

    // J2000 epoch (t1 = 0, simplest case)
    group.bench_function("get_geo_position_j2000", |b| {
        b.iter(|| {
            let jd = black_box(JulianDate::J2000);
            let _p = Moon::get_geo_position::<Kilometer>(jd);
        });
    });

    // Offset epoch (~3.4 years from J2000)
    group.bench_function("get_geo_position_offset", |b| {
        b.iter(|| {
            let jd = black_box(JulianDate::J2000 + Days::new(1234.567));
            let _p = Moon::get_geo_position::<Kilometer>(jd);
        });
    });

    // Far future epoch (~50 years from J2000)
    group.bench_function("get_geo_position_far_future", |b| {
        b.iter(|| {
            let jd = black_box(JulianDate::J2000 + Days::new(18262.5)); // ~50 years
            let _p = Moon::get_geo_position::<Kilometer>(jd);
        });
    });

    // Negative epoch (before J2000)
    group.bench_function("get_geo_position_negative", |b| {
        b.iter(|| {
            let jd = black_box(JulianDate::J2000 - Days::new(3652.5)); // ~10 years before
            let _p = Moon::get_geo_position::<Kilometer>(jd);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_moon_geo_position);
criterion_main!(benches);
