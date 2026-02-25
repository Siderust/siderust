// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compare Sun, Moon, Star Example (prefixed)
//!
//! Run with: `cargo run --example 11_compare_sun_moon_star`

use siderust::bodies::catalog::SIRIUS;
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::time::*;

use qtty::*;

fn main() {
    println!("Generic Altitude Analysis (prefixed)");
}
