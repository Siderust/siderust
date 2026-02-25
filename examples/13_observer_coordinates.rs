// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Observer Coordinates Example (prefixed)
//!
//! Run with: `cargo run --example 13_observer_coordinates`

use qtty::*;
use siderust::coordinates::cartesian::position::EquatorialMeanJ2000;
use siderust::coordinates::centers::{Geocentric, Geodetic};
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical;
use siderust::coordinates::transform::centers::ToTopocentricExt;
use siderust::time::JulianDate;

fn main() {
    println!("Observer coordinates (prefixed)");
}
