// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::{Days, Years};
use siderust::{astro::JulianDate, bodies::solar_system::*, calculus::vsop87::VSOP87};
use std::hint::black_box;

fn bench_vsop87(c: &mut Criterion) {
    let start = JulianDate::J2000 + Years::new(25.0);

    let planet_list: Vec<(&str, Box<dyn VSOP87>)> = vec![
        ("Mercury", Box::new(Mercury)),
        ("Venus", Box::new(Venus)),
        ("Earth", Box::new(Earth)),
        ("Mars", Box::new(Mars)),
        ("Jupiter", Box::new(Jupiter)),
        ("Saturn", Box::new(Saturn)),
        ("Uranus", Box::new(Uranus)),
        ("Neptune", Box::new(Neptune)),
    ];

    for (name, planet) in &planet_list {
        c.bench_function(&format!("vsop87a {}", name), |b| {
            let mut jd = start;
            b.iter(|| {
                jd += Days::new(1.0);
                let _coords = planet.vsop87a(black_box(jd));
            });
        });

        c.bench_function(&format!("vsop87e {}", name), |b| {
            let mut jd = start;
            b.iter(|| {
                jd += Days::new(1.0);
                let _coords = planet.vsop87e(black_box(jd));
            });
        });
    }
}
criterion_group! {
    name = vsop87_benches;
    config = Criterion::default().without_plots();
    targets = bench_vsop87
}
criterion_main!(vsop87_benches);
