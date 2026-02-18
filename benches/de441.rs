// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Micro-benchmarks for DE441 ephemeris body queries.
//!
//! Run with: `cargo bench --bench de441 --features de441`

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::Days;
use siderust::calculus::ephemeris::{De441Ephemeris, Ephemeris};
use siderust::time::JulianDate;
use std::hint::black_box;

fn bench_de441(c: &mut Criterion) {
    let mut group = c.benchmark_group("de441_ephemeris");

    group.bench_function("sun_barycentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::sun_barycentric(black_box(jd));
        });
    });

    group.bench_function("earth_barycentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::earth_barycentric(black_box(jd));
        });
    });

    group.bench_function("earth_heliocentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::earth_heliocentric(black_box(jd));
        });
    });

    group.bench_function("earth_barycentric_velocity", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::earth_barycentric_velocity(black_box(jd));
        });
    });

    group.bench_function("moon_geocentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = De441Ephemeris::moon_geocentric(black_box(jd));
        });
    });

    group.finish();
}

criterion_group!(benches, bench_de441);
criterion_main!(benches);
