// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SPK evaluation and runtime ephemeris benchmarks (synthetic + optional BSP).

use criterion::{criterion_group, criterion_main, Criterion};
use siderust::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::formats::spice::spk::{self, BspSegments, SegmentData};
use siderust::formats::spice::SpkKernelSet;
use siderust::time::JulianDate;
use std::hint::black_box;

const SECONDS_PER_DAY: f64 = siderust::qtty::time::SECONDS_PER_DAY;
const JD_J2000: f64 = tempoch::J2000_JD_TT_DAY.value();

fn synthetic_bsp_segments() -> BspSegments {
    let intlen = 1000.0 * SECONDS_PER_DAY;
    let make = |data_type: i32, rsize: usize, records: Vec<f64>| SegmentData {
        data_type,
        init: 0.0,
        intlen,
        rsize,
        ncoeff: (rsize - 2) / if data_type == 3 { 6 } else { 3 },
        n_records: 1,
        records,
    };
    BspSegments {
        sun: make(
            2,
            8,
            vec![intlen / 2.0, intlen / 2.0, 1.0e8, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        emb: make(
            2,
            8,
            vec![intlen / 2.0, intlen / 2.0, 1.5e8, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        moon: make(
            2,
            8,
            vec![intlen / 2.0, intlen / 2.0, 3.84e5, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        earth: None,
    }
}

fn synthetic_type3_segments() -> BspSegments {
    let intlen = 1000.0 * SECONDS_PER_DAY;
    let mut segs = synthetic_bsp_segments();
    segs.sun = SegmentData {
        data_type: 3,
        init: 0.0,
        intlen,
        rsize: 14,
        ncoeff: 2,
        n_records: 1,
        records: vec![
            intlen / 2.0,
            intlen / 2.0,
            1.0e8,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            2.0,
            0.0,
            3.0,
        ],
    };
    segs
}

fn bench_spk_eval(c: &mut Criterion) {
    let jd = JulianDate::new(JD_J2000 + 500.0);
    let type2 = RuntimeEphemeris::from_segments(synthetic_bsp_segments());
    let type3 = RuntimeEphemeris::from_segments(synthetic_type3_segments());

    let mut group = c.benchmark_group("spk_eval_synthetic");
    group.bench_function("type2_earth_barycentric", |b| {
        b.iter(|| black_box(type2.earth_barycentric(black_box(jd))));
    });
    group.bench_function("type2_earth_barycentric_velocity", |b| {
        b.iter(|| black_box(type2.earth_barycentric_velocity(black_box(jd))));
    });
    group.bench_function("type3_sun_barycentric", |b| {
        b.iter(|| black_box(type3.sun_barycentric(black_box(jd))));
    });
    group.bench_function("type2_moon_geocentric", |b| {
        b.iter(|| black_box(type2.moon_geocentric(black_box(jd))));
    });
    group.finish();
}

fn bench_spk_parse(c: &mut Criterion) {
    let Some(path) = std::env::var("SIDERUST_BSP_PATH").ok() else {
        eprintln!("[spk_perf] SIDERUST_BSP_PATH not set, skipping BSP parse benches");
        return;
    };
    let bytes = std::fs::read(&path).expect("read BSP");
    let mut group = c.benchmark_group("spk_parse");
    group.bench_function("parse_indexed_segments", |b| {
        b.iter(|| black_box(spk::parse_indexed_segments(black_box(&bytes)).unwrap()));
    });
    group.finish();
}

fn bench_runtime_and_kernel_set(c: &mut Criterion) {
    let Some(path) = std::env::var("SIDERUST_BSP_PATH").ok() else {
        eprintln!("[spk_perf] SIDERUST_BSP_PATH not set, skipping runtime BSP benches");
        return;
    };
    let eph = RuntimeEphemeris::from_bsp(&path).expect("load BSP");
    let set = SpkKernelSet::from_paths([&path]).expect("kernel set");
    let mut jd = JulianDate::new(JD_J2000);

    let mut group = c.benchmark_group("runtime_ephemeris_bsp");
    group.bench_function("earth_barycentric", |b| {
        b.iter(|| {
            jd = JulianDate::new(jd.raw().value() + 0.01);
            black_box(eph.earth_barycentric(black_box(jd)));
        });
    });
    group.bench_function("moon_geocentric", |b| {
        b.iter(|| {
            jd = JulianDate::new(jd.raw().value() + 0.01);
            black_box(eph.moon_geocentric(black_box(jd)));
        });
    });
    group.finish();

    let mut group = c.benchmark_group("spk_kernel_set_bsp");
    group.bench_function("geometric_state_earth_moon", |b| {
        b.iter(|| {
            jd = JulianDate::new(jd.raw().value() + 0.01);
            black_box(set.try_geometric_state(301, 399, black_box(jd)).unwrap());
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_spk_eval,
    bench_spk_parse,
    bench_runtime_and_kernel_set
);
criterion_main!(benches);
