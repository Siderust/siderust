// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Coordinate Transformations Example (prefixed)
//!
//! Run with: `cargo run --example 02_coordinate_transformations`

use qtty::*;
use siderust::bodies::solar_system::{Earth, Mars};
use siderust::coordinates::cartesian::position::{
    EclipticMeanJ2000, EquatorialMeanJ2000, GCRS, HCRS, ICRS,
};
use siderust::coordinates::centers::{Geocentric, Heliocentric};
use siderust::coordinates::transform::{Transform, TransformCenter, TransformFrame};
use siderust::time::JulianDate;

fn main() {
    println!("=== Coordinate Transformations Example ===\n");

    let jd = JulianDate::J2000;
    println!("Reference time: J2000.0 (JD {:.1})\n", jd);

    // =========================================================================
    // 1. Frame Transformations (same center)
    // =========================================================================
    println!("1. FRAME TRANSFORMATIONS");
    println!("------------------------");

    // Start with ecliptic coordinates (heliocentric)
    let pos_ecliptic = EclipticMeanJ2000::<Au>::new(1.0, 0.0, 0.0);
    println!("Original (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", pos_ecliptic.x());
    println!("  Y = {:.6} AU", pos_ecliptic.y());
    println!("  Z = {:.6} AU\n", pos_ecliptic.z());

    // Transform to equatorial frame (same heliocentric center)
    let pos_equatorial: EquatorialMeanJ2000<Au, Heliocentric> = pos_ecliptic.to_frame();
    println!("Transformed to EquatorialMeanJ2000 frame:");
    println!("  X = {:.6} AU", pos_equatorial.x());
    println!("  Y = {:.6} AU", pos_equatorial.y());
    println!("  Z = {:.6} AU\n", pos_equatorial.z());

    // Transform to ICRS frame
    let pos_hcrs: HCRS<Au> = pos_equatorial.to_frame();
    println!("Transformed to ICRS frame:");
    println!("  X = {:.6} AU", pos_hcrs.x());
    println!("  Y = {:.6} AU", pos_hcrs.y());
    println!("  Z = {:.6} AU\n", pos_hcrs.z());

    // =========================================================================
    // 2. Center Transformations (same frame)
    // =========================================================================
    println!("2. CENTER TRANSFORMATIONS");
    println!("-------------------------");

    // Get Earth's position (heliocentric ecliptic)
    let earth_helio = Earth::vsop87a(jd);
    println!("Earth (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", earth_helio.x());
    println!("  Y = {:.6} AU", earth_helio.y());
    println!("  Z = {:.6} AU", earth_helio.z());
    println!("  Distance = {:.6} AU\n", earth_helio.distance());

    // Transform to geocentric (Earth becomes origin)
    let earth_geo: EclipticMeanJ2000<Au, Geocentric> = earth_helio.to_center(jd);
    println!("Earth (Geocentric EclipticMeanJ2000) - at origin:");
    println!("  X = {:.10} AU", earth_geo.x());
    println!("  Y = {:.10} AU", earth_geo.y());
    println!("  Z = {:.10} AU\n", earth_geo.z());
    println!(
        "  Distance = {:.10} AU (should be ~0)\n",
        earth_geo.distance()
    );

    // ... (rest of example unchanged) ...
}
