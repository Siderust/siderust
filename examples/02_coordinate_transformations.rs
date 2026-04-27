// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Coordinate Transformations Example (prefixed)
//!
//! Run with: `cargo run --example 02_coordinate_transformations`

use siderust::qtty::*;
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
    println!("  {pos_ecliptic}\n");

    // Transform to equatorial frame (same heliocentric center)
    let pos_equatorial: EquatorialMeanJ2000<Au, Heliocentric> = pos_ecliptic.to_frame();
    println!("Transformed to EquatorialMeanJ2000 frame:");
    println!("  {pos_equatorial}\n");

    // Transform to ICRS frame
    let pos_hcrs: HCRS<Au> = pos_equatorial.to_frame();
    println!("Transformed to ICRS frame:");
    println!("  {pos_hcrs}\n");

    // =========================================================================
    // 2. Center Transformations (same frame)
    // =========================================================================
    println!("2. CENTER TRANSFORMATIONS");
    println!("-------------------------");

    // Get Earth's position (heliocentric ecliptic)
    let earth_helio = Earth::vsop87a(jd);
    println!("Earth (Heliocentric EclipticMeanJ2000):");
    println!("  {earth_helio}");
    println!("  Distance = {:.6}\n", earth_helio.distance());

    // Transform to geocentric (Earth becomes origin)
    let earth_geo: EclipticMeanJ2000<Au, Geocentric> = earth_helio.to_center(jd);
    println!("Earth (Geocentric EclipticMeanJ2000) - at origin:");
    println!("  {earth_geo}");
    println!("  Distance = {:.10} (should be ~0)\n", earth_geo.distance());

    // Get Mars position (heliocentric)
    let mars_helio = Mars::vsop87a(jd);
    println!("Mars (Heliocentric EclipticMeanJ2000):");
    println!("  {mars_helio}");
    println!("  Distance = {:.6}\n", mars_helio.distance());

    // Transform Mars to geocentric
    let mars_geo: EclipticMeanJ2000<Au, Geocentric> = mars_helio.to_center(jd);
    println!("Mars (Geocentric EclipticMeanJ2000) - as seen from Earth:");
    println!("  {mars_geo}");
    println!("  Distance = {:.6}\n", mars_geo.distance());

    // =========================================================================
    // 3. Combined Transformations (center + frame)
    // =========================================================================
    println!("3. COMBINED TRANSFORMATIONS");
    println!("---------------------------");

    // Mars: Heliocentric EclipticMeanJ2000 → Geocentric EquatorialMeanJ2000
    println!("Mars transformation chain:");
    println!("  Start: Heliocentric EclipticMeanJ2000");

    // Method 1: Step by step
    let mars_helio_equ: EquatorialMeanJ2000<Au, Heliocentric> = mars_helio.to_frame();
    println!("  Step 1: Transform frame → Heliocentric EquatorialMeanJ2000");

    let mars_geo_equ: EquatorialMeanJ2000<Au, Geocentric> = mars_helio_equ.to_center(jd);
    println!("  Step 2: Transform center → Geocentric EquatorialMeanJ2000");
    println!("  Result:");
    println!("    {mars_geo_equ}\n");

    // Method 2: Using the Transform trait (does both)
    let mars_geo_equ_direct: EquatorialMeanJ2000<Au, Geocentric> = mars_helio.transform(jd);
    println!("  Or using .transform(jd) directly:");
    println!("    {mars_geo_equ_direct}\n");

    // =========================================================================
    // 4. Barycentric Coordinates
    // =========================================================================
    println!("4. BARYCENTRIC COORDINATES");
    println!("--------------------------");

    // Get Earth in barycentric coordinates
    let earth_bary = Earth::vsop87e(jd);
    println!("Earth (Barycentric EclipticMeanJ2000):");
    println!("  {earth_bary}");
    println!("  Distance from SSB = {:.6}\n", earth_bary.distance());

    // Transform to geocentric
    let earth_geo_from_bary: EclipticMeanJ2000<Au, Geocentric> = earth_bary.to_center(jd);
    println!("Earth (Geocentric, from Barycentric):");
    println!(
        "  Distance = {:.10} (should be ~0)\n",
        earth_geo_from_bary.distance()
    );

    // Transform Mars from barycentric to geocentric
    let mars_bary = Mars::vsop87e(jd);
    let mars_geo_from_bary: EclipticMeanJ2000<Au, Geocentric> = mars_bary.to_center(jd);
    println!("Mars (Geocentric, from Barycentric):");
    println!("  {mars_geo_from_bary}");
    println!("  Distance = {:.6}\n", mars_geo_from_bary.distance());

    // =========================================================================
    // 5. ICRS Frame Transformations
    // =========================================================================
    println!("5. ICRS FRAME TRANSFORMATIONS");
    println!("-----------------------------");

    // Barycentric ICRS (standard for catalogs)
    let star_icrs: ICRS<Au> = ICRS::new(100.0, 50.0, 1000.0);
    println!("Star (Barycentric ICRS):");
    println!("  {star_icrs}\n");

    // Transform to Geocentric ICRS (GCRS)
    let star_gcrs: GCRS<Au> = star_icrs.transform(jd);
    println!("Star (Geocentric ICRS/GCRS):");
    println!("  {star_gcrs}");
    println!("  (Difference is tiny for distant stars)\n");

    // =========================================================================
    // 6. Round-trip Transformation
    // =========================================================================
    println!("6. ROUND-TRIP TRANSFORMATION");
    println!("----------------------------");

    let original = mars_helio;
    println!("Original Mars (Heliocentric EclipticMeanJ2000):");
    println!("  {original}\n");

    // Transform: Helio Ecl → Geo EquatorialMeanJ2000 → Helio Ecl
    let temp: EquatorialMeanJ2000<Au, Geocentric> = original.transform(jd);
    let recovered: EclipticMeanJ2000<Au, Heliocentric> = temp.transform(jd);

    println!("After round-trip transformation:");
    println!("  {recovered}\n");

    let diff_x = (original.x() - recovered.x()).abs();
    let diff_y = (original.y() - recovered.y()).abs();
    let diff_z = (original.z() - recovered.z()).abs();
    println!("Differences (should be tiny):");
    println!("  ΔX = {}", diff_x);
    println!("  ΔY = {}", diff_y);
    println!("  ΔZ = {}\n", diff_z);

    println!("=== Example Complete ===");
}
