// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration tests for `math_core::root_finding` — the astronomy-agnostic
//! Brent & bisection solvers that use typed [`qtty`] quantities.

use qtty::{Day, Quantity, Radian};
use siderust::calculus::math_core::root_finding;
use siderust::time::{ModifiedJulianDate, Period};
use std::cell::Cell;

type Days = Quantity<Day>;
type Mjd = ModifiedJulianDate;
type Radians = Quantity<Radian>;

// =============================================================================
// Brent's method
// =============================================================================

#[test]
fn brent_finds_sine_root_near_pi() {
    let root = root_finding::brent(Days::new(3.0), Days::new(4.0), |t: Days| {
        Radians::new(t.value().sin())
    })
    .expect("should find π");
    assert!(
        (root.value() - std::f64::consts::PI).abs() < 1e-9,
        "expected root near π, got {}",
        root.value()
    );
}

#[test]
fn brent_finds_linear_root() {
    let root = root_finding::brent(Days::new(0.0), Days::new(10.0), |t: Days| {
        Radians::new(t.value() - 5.0)
    })
    .expect("should find 5");
    assert!((root.value() - 5.0).abs() < 1e-10);
}

#[test]
fn brent_handles_step_function() {
    let root = root_finding::brent(Days::new(-2.0), Days::new(2.0), |t: Days| {
        Radians::new(if t < 0.0 { -5.0 } else { 5.0 })
    })
    .expect("step crossing");
    assert!(root.abs() < 1e-6);
}

#[test]
fn brent_returns_none_for_invalid_bracket() {
    assert!(
        root_finding::brent(Days::new(0.0), Days::new(1.0), |_: Days| Radians::new(42.0)).is_none()
    );
}

#[test]
fn brent_returns_endpoint_when_exact() {
    let root = root_finding::brent(Days::new(0.0), Days::new(5.0), |t: Days| {
        Radians::new(t.value() - 5.0)
    })
    .expect("endpoint root");
    assert!((root.value() - 5.0).abs() < 1e-12);
}

#[test]
fn brent_with_values_saves_evaluations() {
    let count = Cell::new(0usize);
    let f = |t: Mjd| -> Radians {
        count.set(count.get() + 1);
        Radians::new(t.value().sin())
    };
    let f_lo = Radians::new((3.0_f64).sin());
    let f_hi = Radians::new((4.0_f64).sin());

    let _ =
        root_finding::brent_with_values(Period::new(Mjd::new(3.0), Mjd::new(4.0)), f_lo, f_hi, f);
    let with_vals = count.get();

    count.set(0);
    let _ = root_finding::brent(Days::new(3.0), Days::new(4.0), |t: Days| {
        count.set(count.get() + 1);
        Radians::new(t.value().sin())
    });
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
        Period::new(Mjd::new(3.0), Mjd::new(4.0)),
        Radians::new((3.0_f64).sin()),
        Radians::new((4.0_f64).sin()),
        |t: Mjd| Radians::new(t.value().sin()),
        Days::new(1e-3),
    )
    .expect("relaxed tolerance");
    assert!((root.value() - std::f64::consts::PI).abs() < 2e-3);
}

#[test]
fn brent_cubic() {
    let root = root_finding::brent(Days::new(1.0), Days::new(2.0), |t: Days| {
        Radians::new(t.value().powi(3) - 2.0)
    })
    .expect("cbrt 2");
    assert!((root.value() - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-9);
}

// =============================================================================
// Bisection
// =============================================================================

#[test]
fn bisection_finds_sine_root() {
    let root = root_finding::bisection(Days::new(3.0), Days::new(4.0), |t: Days| {
        Radians::new(t.value().sin())
    })
    .expect("π via bisection");
    assert!((root.value() - std::f64::consts::PI).abs() < 1e-8);
}

#[test]
fn bisection_returns_none_for_invalid_bracket() {
    assert!(
        root_finding::bisection(Days::new(0.0), Days::new(1.0), |_: Days| Radians::new(42.0))
            .is_none()
    );
}

#[test]
fn bisection_handles_step_function() {
    let root = root_finding::bisection(Days::new(-1.0), Days::new(1.0), |t: Days| {
        Radians::new(if t < 0.0 { -5.0 } else { 5.0 })
    })
    .expect("step crossing via bisection");
    assert!(root.abs() < 1e-6);
}

#[test]
fn bisection_endpoint_root() {
    let root = root_finding::bisection(Days::new(0.0), Days::new(5.0), |t: Days| {
        Radians::new(t.value())
    })
    .expect("root at 0");
    assert!(root.abs() < 1e-12);
}
