// SPDX-License-Identifier: AGPL-3.0-only
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
use siderust::astro::sidereal::{gast_iau2006a, gmst_iau2006};
use siderust::ephemeris::{DynEphemeris, RuntimeEphemeris};
use std::hint::black_box;

fn bench_runtime_ephemeris(c: &mut Criterion) {
    let Some(path) = std::env::var_os("SIDERUST_BSP_PATH") else {
        eprintln!("[de441 bench] SIDERUST_BSP_PATH not set, skipping RuntimeEphemeris benchmarks");
        return;
    };

    let eph = RuntimeEphemeris::from_bsp(&path)
        .unwrap_or_else(|err| panic!("failed to load BSP {path:?}: {err}"));
    let jd = siderust::time::J2000;
    let mut group = c.benchmark_group("runtime_ephemeris");

    group.bench_function("sun_barycentric", |b| {
        eph.try_sun_barycentric(jd)
            .expect("BSP must cover J2000 TT for sun_barycentric");
        b.iter(|| eph.sun_barycentric(black_box(jd)));
    });

    group.bench_function("earth_barycentric", |b| {
        eph.try_earth_barycentric(jd)
            .expect("BSP must cover J2000 TT for earth_barycentric");
        b.iter(|| eph.earth_barycentric(black_box(jd)));
    });

    group.bench_function("earth_heliocentric", |b| {
        eph.try_earth_heliocentric(jd)
            .expect("BSP must cover J2000 TT for earth_heliocentric");
        b.iter(|| eph.earth_heliocentric(black_box(jd)));
    });

    group.bench_function("earth_barycentric_velocity", |b| {
        eph.try_earth_barycentric_velocity(jd)
            .expect("BSP must cover J2000 TT for earth_barycentric_velocity");
        b.iter(|| eph.earth_barycentric_velocity(black_box(jd)));
    });

    group.bench_function("moon_geocentric", |b| {
        eph.try_moon_geocentric(jd)
            .expect("BSP must cover J2000 TT for moon_geocentric");
        b.iter(|| eph.moon_geocentric(black_box(jd)));
    });

    group.bench_function("gmst_iau2006", |b| {
        let jd = siderust::time::J2000;
        let ut1 = jd
            .to_with::<tempoch::UT1>(&tempoch::TimeContext::new())
            .unwrap();
        b.iter(|| black_box(gmst_iau2006(black_box(ut1), black_box(jd))));
    });

    group.bench_function("gast_iau2006a", |b| {
        let jd = siderust::time::J2000;
        let ut1 = jd
            .to_with::<tempoch::UT1>(&tempoch::TimeContext::new())
            .unwrap();
        b.iter(|| black_box(gast_iau2006a(black_box(ut1), black_box(jd))));
    });

    group.finish();
}

criterion_group!(benches, bench_runtime_ephemeris);
criterion_main!(benches);
