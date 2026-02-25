// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! All Center Conversions Example (prefixed)
//!
//! Run with: `cargo run --example 04_all_center_conversions`

use qtty::*;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::{Barycentric, Geocentric, Heliocentric, ReferenceCenter};
use siderust::coordinates::frames::EclipticMeanJ2000;
use siderust::coordinates::transform::{CenterShiftProvider, PositionAstroExt};
use siderust::time::JulianDate;

type F = EclipticMeanJ2000;
type U = AstronomicalUnit;

fn main() {
    let jd = JulianDate::new(2_460_000.5);
    println!("Center conversion demo at JD(TT) = {:.1}", jd.value());
}
