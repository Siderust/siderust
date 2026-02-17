// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Brent Root Finding Example
//!
//! Demonstrates:
//! - `math_core::root_finding::brent` on typed `qtty` quantities
//! - Unit-safe root finding on `f(x) = 0` within a bracket `[lo, hi]`

use qtty::{Quantity, Radians, Unitless};
use siderust::calculus::math_core::root_finding::brent;

fn main() {
    println!("=== Siderust Brent Root Finding Example ===\n");

    // Solve sin(x) - 0.5 = 0 on [0, 2] radians.
    let lo = Radians::new(0.0);
    let hi = Radians::new(2.0);

    let root = brent(lo, hi, |x: Radians| Quantity::<Unitless>::new(x.sin() - 0.5))
        .expect("bracket does not contain a root");

    let residual = (root.sin() - 0.5).abs();

    println!("Root for sin(x) = 0.5 in [0, 2] rad:");
    println!("  x ≈ {}", root);
    println!("  residual = {:.3e}", residual);
}

