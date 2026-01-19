// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::{Au, Degrees};
use siderust::coordinates::{cartesian, centers::Barycentric, frames::ICRS, spherical};
use std::{hint::black_box, time::Duration};

fn bench_cartesian_spherical_converters(c: &mut Criterion) {
    // Test data: one Cartesian and one Spherical point in the ICRS/Barycentric frame
    let icrs_cartesian = cartesian::Position::<Barycentric, ICRS, Au>::new(10.0, 20.0, 30.0);
    let icrs_spherical = spherical::Position::<Barycentric, ICRS, Au>::new(
        Degrees::new(10.0),
        Degrees::new(20.0),
        30.0,
    );

    // Create a benchmark group so we can tweak settings just for these benchmarks
    let mut group = c.benchmark_group("converters");

    // Collect more samples and run a bit longer to reduce statistical noise
    group
        .sample_size(1_000) // number of samples to collect
        .measurement_time(Duration::from_secs(3)) // time spent measuring
        .warm_up_time(Duration::from_secs(1)); // time to let caches/freq settle

    // Cartesian → Spherical
    group.bench_function("to_spherical (ICRS)", |b| {
        b.iter(|| {
            // black_box prevents the compiler from optimizing the value away
            let cart = black_box(&icrs_cartesian);
            let res = cart.to_spherical();
            black_box(res); // keep the result “alive”
        });
    });

    // Spherical → Cartesian
    group.bench_function("to_cartesian (ICRS)", |b| {
        b.iter(|| {
            let sph = black_box(&icrs_spherical);
            let res = sph.to_cartesian();
            black_box(res);
        });
    });

    group.finish();
}

// Disable plot generation (faster CI runs); keep default configuration otherwise
criterion_group! {
    name = converter_benches;
    config = Criterion::default().without_plots();
    targets = bench_cartesian_spherical_converters
}
criterion_main!(converter_benches);
