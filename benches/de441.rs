// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Micro-benchmarks for runtime DE4xx ephemeris body queries.
//!
//! Requires a BSP file; skips gracefully when none is present.
//!
//! Run with:
//! ```bash
//! SIDERUST_BSP_PATH=/path/to/de440.bsp cargo bench --bench de441
//! ```

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::Days;
use siderust::calculus::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::time::JulianDate;
use std::hint::black_box;

fn bench_runtime_ephemeris(c: &mut Criterion) {
    let path = match std::env::var("SIDERUST_BSP_PATH") {
        Ok(p) => p,
        Err(_) => {
            eprintln!(
                "[de441 bench] SIDERUST_BSP_PATH not set — skipping RuntimeEphemeris benchmarks"
            );
            return;
        }
    };

    let eph = match RuntimeEphemeris::from_bsp(&path) {
        Ok(e) => e,
        Err(err) => {
            eprintln!("[de441 bench] Failed to load '{path}': {err} — skipping");
            return;
        }
    };

    let mut group = c.benchmark_group("runtime_ephemeris");

    group.bench_function("sun_barycentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = eph.sun_barycentric(black_box(jd));
        });
    });

    group.bench_function("earth_barycentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = eph.earth_barycentric(black_box(jd));
        });
    });

    group.bench_function("earth_heliocentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = eph.earth_heliocentric(black_box(jd));
        });
    });

    group.bench_function("earth_barycentric_velocity", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = eph.earth_barycentric_velocity(black_box(jd));
        });
    });

    group.bench_function("moon_geocentric", |b| {
        let mut jd = JulianDate::J2000;
        b.iter(|| {
            jd += Days::new(1.0);
            let _ = eph.moon_geocentric(black_box(jd));
        });
    });

    group.finish();
}

criterion_group!(benches, bench_runtime_ephemeris);
criterion_main!(benches);
