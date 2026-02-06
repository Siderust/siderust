// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for `math_core::root_finding` — the astronomy-agnostic
//! Brent & bisection solvers that replaced the old JulianDate-typed module.

use siderust::calculus::math_core::root_finding;
use std::cell::Cell;

// =============================================================================
// Brent's method
// =============================================================================

#[test]
fn brent_finds_sine_root_near_pi() {
    let root = root_finding::brent(3.0, 4.0, |t| t.sin()).expect("should find π");
    assert!(
        (root - std::f64::consts::PI).abs() < 1e-9,
        "expected root near π, got {}",
        root
    );
}

#[test]
fn brent_finds_linear_root() {
    let root = root_finding::brent(0.0, 10.0, |t| t - 5.0).expect("should find 5");
    assert!((root - 5.0).abs() < 1e-10);
}

#[test]
fn brent_handles_step_function() {
    let root = root_finding::brent(-2.0, 2.0, |t| if t < 0.0 { -5.0 } else { 5.0 })
        .expect("step crossing");
    assert!(root.abs() < 1e-6);
}

#[test]
fn brent_returns_none_for_invalid_bracket() {
    assert!(root_finding::brent(0.0, 1.0, |_| 42.0).is_none());
}

#[test]
fn brent_returns_endpoint_when_exact() {
    let root = root_finding::brent(0.0, 5.0, |t| t - 5.0).expect("endpoint root");
    assert!((root - 5.0).abs() < 1e-12);
}

#[test]
fn brent_with_values_saves_evaluations() {
    let count = Cell::new(0usize);
    let f = |t: f64| {
        count.set(count.get() + 1);
        t.sin()
    };
    let f_lo = (3.0_f64).sin();
    let f_hi = (4.0_f64).sin();

    let _ = root_finding::brent_with_values(3.0, 4.0, f_lo, f_hi, &f);
    let with_vals = count.get();

    count.set(0);
    let _ = root_finding::brent(3.0, 4.0, &f);
    let without = count.get();

    // brent_with_values saves the 2 endpoint evaluations
    assert!(
        with_vals <= without,
        "with_values used {} calls, plain brent used {}",
        with_vals,
        without
    );
}

#[test]
fn brent_tol_respects_relaxed_tolerance() {
    let root = root_finding::brent_tol(
        3.0,
        4.0,
        (3.0_f64).sin(),
        (4.0_f64).sin(),
        |t| t.sin(),
        1e-3,
    )
    .expect("relaxed tolerance");
    assert!((root - std::f64::consts::PI).abs() < 2e-3);
}

#[test]
fn brent_cubic() {
    let root = root_finding::brent(1.0, 2.0, |t| t.powi(3) - 2.0).expect("cbrt 2");
    assert!((root - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-9);
}

// =============================================================================
// Bisection
// =============================================================================

#[test]
fn bisection_finds_sine_root() {
    let root = root_finding::bisection(3.0, 4.0, |t| t.sin()).expect("π via bisection");
    assert!((root - std::f64::consts::PI).abs() < 1e-8);
}

#[test]
fn bisection_returns_none_for_invalid_bracket() {
    assert!(root_finding::bisection(0.0, 1.0, |_| 42.0).is_none());
}

#[test]
fn bisection_handles_step_function() {
    let root = root_finding::bisection(-1.0, 1.0, |t| if t < 0.0 { -5.0 } else { 5.0 })
        .expect("step crossing via bisection");
    assert!(root.abs() < 1e-6);
}

#[test]
fn bisection_endpoint_root() {
    let root = root_finding::bisection(0.0, 5.0, |t| t).expect("root at 0");
    assert!(root.abs() < 1e-12);
}
