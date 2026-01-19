// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use siderust::astro::JulianDate;
use siderust::calculus::root_finding::{find_crossing, refine_root_bisection, refine_root_newton};
use std::cell::Cell;

#[test]
fn finds_sine_root_near_pi() {
    let lo = JulianDate::new(3.0);
    let hi = JulianDate::new(4.0);
    let scalar = |jd: JulianDate| jd.value().sin();

    let result = find_crossing(lo, hi, &scalar, 0.0).expect("should find pi");
    assert!((result.value() - std::f64::consts::PI).abs() < 1e-9);
}

#[test]
fn step_function_crossing_returns_zero() {
    let lo = JulianDate::new(-2.0);
    let hi = JulianDate::new(2.0);
    let scalar = |jd: JulianDate| {
        if jd.value() < 0.0 {
            -5.0
        } else {
            5.0
        }
    };

    let result = find_crossing(lo, hi, &scalar, 0.0).expect("step crossing should succeed");
    assert!(result.value().abs() < 1e-6);
}

#[test]
fn bisection_endpoint_match_returns_hi() {
    let threshold = 10_000.0;
    let lo = JulianDate::new(threshold - 10.0);
    let hi = JulianDate::new(threshold);
    let scalar = |jd: JulianDate| jd.value();

    let result = refine_root_bisection(lo, hi, scalar, threshold);
    assert_eq!(result, Some(hi));
}

#[test]
fn newton_hits_max_iterations_with_stateful_scalar() {
    let calls = Cell::new(0usize);
    let scalar = |_jd: JulianDate| {
        let count = calls.get() + 1;
        calls.set(count);
        count as f64
    };

    let guess = JulianDate::J2000;
    assert!(refine_root_newton(guess, scalar, 0.0).is_none());
}
