// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Coordinate Transformations Example
//!
//! This example demonstrates transformations between different coordinate systems:
//! - Frame transformations (EclipticMeanJ2000 ↔ EquatorialMeanJ2000 ↔ ICRS)
//! - Center transformations (Heliocentric ↔ Geocentric ↔ Barycentric)
//! - Combined transformations
//! - Time-dependent transformations

use qtty::*;
use siderust::bodies::solar_system::{Earth, Mars};
use siderust::coordinates::cartesian::position::{EclipticMeanJ2000, EquatorialMeanJ2000, GCRS, HCRS, ICRS};
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
    let earth_helio = *Earth::vsop87a(jd).get_position();
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
    println!("  Z = {:.10} AU", earth_geo.z());
    println!(
        "  Distance = {:.10} AU (should be ~0)\n",
        earth_geo.distance()
    );

    // Get Mars position (heliocentric)
    let mars_helio = *Mars::vsop87a(jd).get_position();
    println!("Mars (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", mars_helio.x());
    println!("  Y = {:.6} AU", mars_helio.y());
    println!("  Z = {:.6} AU", mars_helio.z());
    println!("  Distance = {:.6} AU\n", mars_helio.distance());

    // Transform Mars to geocentric
    let mars_geo: EclipticMeanJ2000<Au, Geocentric> = mars_helio.to_center(jd);
    println!("Mars (Geocentric EclipticMeanJ2000) - as seen from Earth:");
    println!("  X = {:.6} AU", mars_geo.x());
    println!("  Y = {:.6} AU", mars_geo.y());
    println!("  Z = {:.6} AU", mars_geo.z());
    println!("  Distance = {:.6} AU\n", mars_geo.distance());

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
    println!("    X = {:.6} AU", mars_geo_equ.x());
    println!("    Y = {:.6} AU", mars_geo_equ.y());
    println!("    Z = {:.6} AU\n", mars_geo_equ.z());

    // Method 2: Using the Transform trait (does both)
    let mars_geo_equ_direct: EquatorialMeanJ2000<Au, Geocentric> = mars_helio.transform(jd);
    println!("  Or using .transform(jd) directly:");
    println!("    X = {:.6} AU", mars_geo_equ_direct.x());
    println!("    Y = {:.6} AU", mars_geo_equ_direct.y());
    println!("    Z = {:.6} AU\n", mars_geo_equ_direct.z());

    // =========================================================================
    // 4. Barycentric Coordinates
    // =========================================================================
    println!("4. BARYCENTRIC COORDINATES");
    println!("--------------------------");

    // Get Earth in barycentric coordinates
    let earth_bary = *Earth::vsop87e(jd).get_position();
    println!("Earth (Barycentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", earth_bary.x());
    println!("  Y = {:.6} AU", earth_bary.y());
    println!("  Z = {:.6} AU", earth_bary.z());
    println!("  Distance from SSB = {:.6} AU\n", earth_bary.distance());

    // Transform to geocentric
    let earth_geo_from_bary: EclipticMeanJ2000<Au, Geocentric> = earth_bary.to_center(jd);
    println!("Earth (Geocentric, from Barycentric):");
    println!(
        "  Distance = {:.10} AU (should be ~0)\n",
        earth_geo_from_bary.distance()
    );

    // Transform Mars from barycentric to geocentric
    let mars_bary = *Mars::vsop87e(jd).get_position();
    let mars_geo_from_bary: EclipticMeanJ2000<Au, Geocentric> = mars_bary.to_center(jd);
    println!("Mars (Geocentric, from Barycentric):");
    println!("  X = {:.6} AU", mars_geo_from_bary.x());
    println!("  Y = {:.6} AU", mars_geo_from_bary.y());
    println!("  Z = {:.6} AU", mars_geo_from_bary.z());
    println!("  Distance = {:.6} AU\n", mars_geo_from_bary.distance());

    // =========================================================================
    // 5. ICRS Frame Transformations
    // =========================================================================
    println!("5. ICRS FRAME TRANSFORMATIONS");
    println!("-----------------------------");

    // Barycentric ICRS (standard for catalogs)
    let star_icrs: ICRS<Au> = ICRS::new(100.0, 50.0, 1000.0);
    println!("Star (Barycentric ICRS):");
    println!("  X = {:.3} AU", star_icrs.x());
    println!("  Y = {:.3} AU", star_icrs.y());
    println!("  Z = {:.3} AU\n", star_icrs.z());

    // Transform to Geocentric ICRS (GCRS)
    let star_gcrs: GCRS<Au> = star_icrs.transform(jd);
    println!("Star (Geocentric ICRS/GCRS):");
    println!("  X = {:.3} AU", star_gcrs.x());
    println!("  Y = {:.3} AU", star_gcrs.y());
    println!("  Z = {:.3} AU", star_gcrs.z());
    println!("  (Difference is tiny for distant stars)\n");

    // =========================================================================
    // 6. Round-trip Transformation
    // =========================================================================
    println!("6. ROUND-TRIP TRANSFORMATION");
    println!("----------------------------");

    let original = mars_helio;
    println!("Original Mars (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.10} AU", original.x());
    println!("  Y = {:.10} AU", original.y());
    println!("  Z = {:.10} AU\n", original.z());

    // Transform: Helio Ecl → Geo EquatorialMeanJ2000 → Helio Ecl
    let temp: EquatorialMeanJ2000<Au, Geocentric> = original.transform(jd);
    let recovered: EclipticMeanJ2000<Au, Heliocentric> = temp.transform(jd);

    println!("After round-trip transformation:");
    println!("  X = {:.10} AU", recovered.x());
    println!("  Y = {:.10} AU", recovered.y());
    println!("  Z = {:.10} AU\n", recovered.z());

    let diff_x = (original.x() - recovered.x()).abs();
    let diff_y = (original.y() - recovered.y()).abs();
    let diff_z = (original.z() - recovered.z()).abs();
    println!("Differences (should be tiny):");
    println!("  ΔX = {}", diff_x);
    println!("  ΔY = {}", diff_y);
    println!("  ΔZ = {}\n", diff_z);

    println!("=== Example Complete ===");
}
