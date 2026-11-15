// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Micro-benchmarks for ELP2000 evaluation.
//!
//! Run with: `cargo bench --bench elp2000`

use criterion::{criterion_group, criterion_main, Criterion};
use siderust::bodies::solar_system::Moon;
use siderust::qtty::{Days, Kilometer};
use std::hint::black_box;

fn bench_moon_geo_position(c: &mut Criterion) {
    let mut group = c.benchmark_group("elp2000_moon_geo");

    // J2000 epoch (t1 = 0, simplest case)
    group.bench_function("get_geo_position_j2000", |b| {
        let jd = siderust::time::J2000;
        b.iter(|| Moon::get_geo_position::<Kilometer>(black_box(jd)));
    });

    // Offset epoch (~3.4 years from J2000)
    group.bench_function("get_geo_position_offset", |b| {
        let jd = siderust::time::J2000 + Days::new(1234.567);
        b.iter(|| Moon::get_geo_position::<Kilometer>(black_box(jd)));
    });

    // Far future epoch (~50 years from J2000)
    group.bench_function("get_geo_position_far_future", |b| {
        let jd = siderust::time::J2000 + Days::new(18262.5);
        b.iter(|| Moon::get_geo_position::<Kilometer>(black_box(jd)));
    });

    // Negative epoch (before J2000)
    group.bench_function("get_geo_position_negative", |b| {
        let jd = siderust::time::J2000 - Days::new(3652.5);
        b.iter(|| Moon::get_geo_position::<Kilometer>(black_box(jd)));
    });

    group.finish();
}

criterion_group!(benches, bench_moon_geo_position);
criterion_main!(benches);
