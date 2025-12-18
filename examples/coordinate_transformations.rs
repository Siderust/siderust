//! Coordinate Transformations Example
//!
//! This example demonstrates transformations between different coordinate systems:
//! - Frame transformations (Ecliptic ↔ Equatorial ↔ ICRS)
//! - Center transformations (Heliocentric ↔ Geocentric ↔ Barycentric)
//! - Combined transformations
//! - Time-dependent transformations

use siderust::astro::JulianDate;
use siderust::bodies::solar_system::{Earth, Mars};
use siderust::coordinates::cartesian::position::{Ecliptic, Equatorial, GCRS, HCRS, ICRS};
use siderust::coordinates::centers::{Geocentric, Heliocentric};
use siderust::coordinates::transform::{Transform, TransformFrame, TransformCenter};
use qtty::*;

fn main() {
    println!("=== Coordinate Transformations Example ===\n");

    let jd = JulianDate::J2000;
    println!("Reference time: J2000.0 (JD {})\n", jd.value());

    // =========================================================================
    // 1. Frame Transformations (same center)
    // =========================================================================
    println!("1. FRAME TRANSFORMATIONS");
    println!("------------------------");

    // Start with ecliptic coordinates (heliocentric)
    let pos_ecliptic = Ecliptic::<Au>::new(1.0, 0.0, 0.0);
    println!("Original (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", pos_ecliptic.x().value());
    println!("  Y = {:.6} AU", pos_ecliptic.y().value());
    println!("  Z = {:.6} AU\n", pos_ecliptic.z().value());

    // Transform to equatorial frame (same heliocentric center)
    let pos_equatorial: Equatorial<Au, Heliocentric> = pos_ecliptic.to_frame();
    println!("Transformed to Equatorial frame:");
    println!("  X = {:.6} AU", pos_equatorial.x().value());
    println!("  Y = {:.6} AU", pos_equatorial.y().value());
    println!("  Z = {:.6} AU\n", pos_equatorial.z().value());

    // Transform to ICRS frame
    let pos_hcrs: HCRS<Au> = pos_equatorial.to_frame();
    println!("Transformed to ICRS frame:");
    println!("  X = {:.6} AU", pos_hcrs.x().value());
    println!("  Y = {:.6} AU", pos_hcrs.y().value());
    println!("  Z = {:.6} AU\n", pos_hcrs.z().value());

    // =========================================================================
    // 2. Center Transformations (same frame)
    // =========================================================================
    println!("2. CENTER TRANSFORMATIONS");
    println!("-------------------------");

    // Get Earth's position (heliocentric ecliptic)
    let earth_helio = *Earth::vsop87a(jd).get_position();
    println!("Earth (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", earth_helio.x().value());
    println!("  Y = {:.6} AU", earth_helio.y().value());
    println!("  Z = {:.6} AU", earth_helio.z().value());
    println!("  Distance = {:.6} AU\n", earth_helio.distance().value());

    // Transform to geocentric (Earth becomes origin)
    let earth_geo: Ecliptic<Au, Geocentric> = earth_helio.to_center(jd);
    println!("Earth (Geocentric Ecliptic) - at origin:");
    println!("  X = {:.10} AU", earth_geo.x().value());
    println!("  Y = {:.10} AU", earth_geo.y().value());
    println!("  Z = {:.10} AU", earth_geo.z().value());
    println!("  Distance = {:.10} AU (should be ~0)\n", earth_geo.distance().value());

    // Get Mars position (heliocentric)
    let mars_helio = *Mars::vsop87a(jd).get_position();
    println!("Mars (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", mars_helio.x().value());
    println!("  Y = {:.6} AU", mars_helio.y().value());
    println!("  Z = {:.6} AU", mars_helio.z().value());
    println!("  Distance = {:.6} AU\n", mars_helio.distance().value());

    // Transform Mars to geocentric
    let mars_geo: Ecliptic<Au, Geocentric> = mars_helio.to_center(jd);
    println!("Mars (Geocentric Ecliptic) - as seen from Earth:");
    println!("  X = {:.6} AU", mars_geo.x().value());
    println!("  Y = {:.6} AU", mars_geo.y().value());
    println!("  Z = {:.6} AU", mars_geo.z().value());
    println!("  Distance = {:.6} AU\n", mars_geo.distance().value());

    // =========================================================================
    // 3. Combined Transformations (center + frame)
    // =========================================================================
    println!("3. COMBINED TRANSFORMATIONS");
    println!("---------------------------");

    // Mars: Heliocentric Ecliptic → Geocentric Equatorial
    println!("Mars transformation chain:");
    println!("  Start: Heliocentric Ecliptic");
    
    // Method 1: Step by step
    let mars_helio_equ: Equatorial<Au, Heliocentric> = mars_helio.to_frame();
    println!("  Step 1: Transform frame → Heliocentric Equatorial");
    
    let mars_geo_equ: Equatorial<Au, Geocentric> = mars_helio_equ.to_center(jd);
    println!("  Step 2: Transform center → Geocentric Equatorial");
    println!("  Result:");
    println!("    X = {:.6} AU", mars_geo_equ.x().value());
    println!("    Y = {:.6} AU", mars_geo_equ.y().value());
    println!("    Z = {:.6} AU\n", mars_geo_equ.z().value());

    // Method 2: Using the Transform trait (does both)
    let mars_geo_equ_direct: Equatorial<Au, Geocentric> = mars_helio.transform(jd);
    println!("  Or using .transform(jd) directly:");
    println!("    X = {:.6} AU", mars_geo_equ_direct.x().value());
    println!("    Y = {:.6} AU", mars_geo_equ_direct.y().value());
    println!("    Z = {:.6} AU\n", mars_geo_equ_direct.z().value());

    // =========================================================================
    // 4. Barycentric Coordinates
    // =========================================================================
    println!("4. BARYCENTRIC COORDINATES");
    println!("--------------------------");

    // Get Earth in barycentric coordinates
    let earth_bary = *Earth::vsop87e(jd).get_position();
    println!("Earth (Barycentric Ecliptic):");
    println!("  X = {:.6} AU", earth_bary.x().value());
    println!("  Y = {:.6} AU", earth_bary.y().value());
    println!("  Z = {:.6} AU", earth_bary.z().value());
    println!("  Distance from SSB = {:.6} AU\n", earth_bary.distance().value());

    // Transform to geocentric
    let earth_geo_from_bary: Ecliptic<Au, Geocentric> = earth_bary.to_center(jd);
    println!("Earth (Geocentric, from Barycentric):");
    println!("  Distance = {:.10} AU (should be ~0)\n", earth_geo_from_bary.distance().value());

    // Transform Mars from barycentric to geocentric
    let mars_bary = *Mars::vsop87e(jd).get_position();
    let mars_geo_from_bary: Ecliptic<Au, Geocentric> = mars_bary.to_center(jd);
    println!("Mars (Geocentric, from Barycentric):");
    println!("  X = {:.6} AU", mars_geo_from_bary.x().value());
    println!("  Y = {:.6} AU", mars_geo_from_bary.y().value());
    println!("  Z = {:.6} AU", mars_geo_from_bary.z().value());
    println!("  Distance = {:.6} AU\n", mars_geo_from_bary.distance().value());

    // =========================================================================
    // 5. ICRS Frame Transformations
    // =========================================================================
    println!("5. ICRS FRAME TRANSFORMATIONS");
    println!("-----------------------------");

    // Barycentric ICRS (standard for catalogs)
    let star_icrs: ICRS<Au> = ICRS::new(100.0, 50.0, 1000.0);
    println!("Star (Barycentric ICRS):");
    println!("  X = {:.3} AU", star_icrs.x().value());
    println!("  Y = {:.3} AU", star_icrs.y().value());
    println!("  Z = {:.3} AU\n", star_icrs.z().value());

    // Transform to Geocentric ICRS (GCRS)
    let star_gcrs: GCRS<Au> = star_icrs.transform(jd);
    println!("Star (Geocentric ICRS/GCRS):");
    println!("  X = {:.3} AU", star_gcrs.x().value());
    println!("  Y = {:.3} AU", star_gcrs.y().value());
    println!("  Z = {:.3} AU", star_gcrs.z().value());
    println!("  (Difference is tiny for distant stars)\n");

    // =========================================================================
    // 6. Round-trip Transformation
    // =========================================================================
    println!("6. ROUND-TRIP TRANSFORMATION");
    println!("----------------------------");

    let original = mars_helio;
    println!("Original Mars (Heliocentric Ecliptic):");
    println!("  X = {:.10} AU", original.x().value());
    println!("  Y = {:.10} AU", original.y().value());
    println!("  Z = {:.10} AU\n", original.z().value());

    // Transform: Helio Ecl → Geo Equ → Helio Ecl
    let temp: Equatorial<Au, Geocentric> = original.transform(jd);
    let recovered: Ecliptic<Au, Heliocentric> = temp.transform(jd);
    
    println!("After round-trip transformation:");
    println!("  X = {:.10} AU", recovered.x().value());
    println!("  Y = {:.10} AU", recovered.y().value());
    println!("  Z = {:.10} AU\n", recovered.z().value());

    let diff_x = (original.x().value() - recovered.x().value()).abs();
    let diff_y = (original.y().value() - recovered.y().value()).abs();
    let diff_z = (original.z().value() - recovered.z().value()).abs();
    println!("Differences (should be tiny):");
    println!("  ΔX = {:.2e} AU", diff_x);
    println!("  ΔY = {:.2e} AU", diff_y);
    println!("  ΔZ = {:.2e} AU\n", diff_z);

    println!("=== Example Complete ===");
}
