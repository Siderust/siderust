// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Sidereal time benchmarks (GMST / GAST).

use criterion::{criterion_group, criterion_main, Criterion};
use siderust::astro::sidereal::{gast_iau2006, gast_iau2006a, gmst_iau2006};
use siderust::time::JulianDate;
use std::hint::black_box;

fn bench_sidereal(c: &mut Criterion) {
    let jd = JulianDate::new(2_460_000.5);
    let mut group = c.benchmark_group("sidereal");

    group.bench_function("gmst_iau2006", |b| {
        b.iter(|| black_box(gmst_iau2006(black_box(jd), jd)));
    });
    group.bench_function("gast_iau2006a", |b| {
        b.iter(|| black_box(gast_iau2006a(black_box(jd), jd)));
    });
    group.bench_function("gast_iau2006_low_level", |b| {
        let nut = siderust::astro::nutation::nutation_iau2000b(jd);
        b.iter(|| {
            black_box(gast_iau2006(
                black_box(jd),
                jd,
                nut.dpsi,
                nut.mean_obliquity,
            ))
        });
    });
    group.finish();
}

criterion_group!(benches, bench_sidereal);
criterion_main!(benches);
