// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! JPL Precise Ephemeris Example (prefixed)
//!
//! Run with: `cargo run --example 16_jpl_precise_ephemeris --features de440`

#[cfg(not(any(feature = "de440", feature = "de441")))]
fn main() {
    eprintln!("ERROR: This example requires the 'de440' and/or 'de441' feature.");
    std::process::exit(1);
}

#[cfg(any(feature = "de440", feature = "de441"))]
fn main() {
    println!("JPL ephemeris example (prefixed)");
}
