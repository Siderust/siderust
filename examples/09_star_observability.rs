// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Star Observability Planner (prefixed)
//!
//! Run with: `cargo run --example 09_star_observability`

use siderust::bodies::catalog::{ALTAIR, BETELGEUSE, POLARIS, RIGEL, SIRIUS, VEGA};
use siderust::bodies::solar_system::Sun;
use siderust::bodies::Star;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

fn main() {
    println!("Star observability (prefixed)");
}
