// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Targets + Proper Motion Example (prefixed)
//!
//! Run with: `cargo run --example 15_targets_proper_motion`

use qtty::*;
use siderust::astro::proper_motion::{set_proper_motion_since_j2000, ProperMotion};
use siderust::coordinates::spherical::position;
use siderust::targets::CoordinateWithPM;
use siderust::time::JulianDate;

fn main() {
    println!("Targets + proper motion (prefixed)");
}
