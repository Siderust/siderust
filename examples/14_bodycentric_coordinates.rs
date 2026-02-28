// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Bodycentric Coordinates Example (prefixed)
//!
//! Run with: `cargo run --example 14_bodycentric_coordinates`

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::bodies::solar_system::{Earth, Mars, Venus};
use siderust::coordinates::cartesian::position::{EclipticMeanJ2000, Position};
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::centers::{Bodycentric, BodycentricParams, Geocentric, Heliocentric};
use siderust::coordinates::frames;
use siderust::coordinates::transform::TransformCenter;
use siderust::time::JulianDate;

fn main() {
    println!("Bodycentric coordinates (prefixed)");
}
