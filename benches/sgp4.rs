// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Comparable public API initialization and propagation benchmarks.
use criterion::{criterion_group, criterion_main, Criterion};
use siderust::astro::sgp4::TlePropagator;
use siderust::formats::tle::parse_tle;
use siderust::qtty::Minutes;
use std::hint::black_box;

fn benchmark(c: &mut Criterion) {
    let lines: Vec<_> = include_str!("../tests/fixtures/sgp4/vallado.txt")
        .lines()
        .collect();
    let mut group = c.benchmark_group("sgp4");
    group
        .sample_size(30)
        .warm_up_time(std::time::Duration::from_secs(1))
        .measurement_time(std::time::Duration::from_secs(2));
    // Verification satellites: drag-dominated LEO, synchronous GEO, half-day resonance.
    for (name, id) in [("leo", "06251"), ("geo", "24208"), ("half_day", "08195")] {
        let index = lines
            .iter()
            .position(|l| l.starts_with(&format!("1 {id}")))
            .unwrap();
        let tle = parse_tle(lines[index], lines[index + 1]).unwrap();
        group.bench_function(format!("{name}/initialize"), |b| {
            b.iter(|| black_box(TlePropagator::from_tle(black_box(&tle)).unwrap()))
        });
        let propagator = TlePropagator::from_tle(&tle).unwrap();
        group.bench_function(format!("{name}/propagate"), |b| {
            b.iter(|| {
                black_box(
                    propagator
                        .propagate_minutes(black_box(Minutes::new(1440.25)))
                        .unwrap(),
                )
            })
        });
        group.bench_function(format!("{name}/repeated_100"), |b| {
            b.iter(|| {
                for i in 0..100 {
                    black_box(
                        propagator
                            .propagate_minutes(Minutes::new(i as f64 * 60.25))
                            .unwrap(),
                    );
                }
            })
        });
    }
    group.finish();
}
criterion_group!(benches, benchmark);
criterion_main!(benches);
