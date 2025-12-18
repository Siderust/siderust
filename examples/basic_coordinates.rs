//! Basic Coordinates Example
//!
//! This example demonstrates the fundamental concepts of the siderust coordinate system:
//! - Creating positions and directions
//! - Different reference frames and centers
//! - Cartesian vs Spherical representations
//! - Type safety and units

use siderust::coordinates::cartesian;
use siderust::coordinates::spherical;
use siderust::coordinates::centers::{self, ReferenceCenter, ObserverSite};
use siderust::coordinates::frames::{self, ReferenceFrame};
use qtty::*;

fn main() {
    println!("=== Siderust Basic Coordinates Example ===\n");

    // =========================================================================
    // 1. Cartesian Coordinates
    // =========================================================================
    println!("1. CARTESIAN COORDINATES");
    println!("------------------------");

    // Create a heliocentric ecliptic position (1 AU along X-axis)
    let earth_position = cartesian::position::Ecliptic::<AstronomicalUnit>::new(
        1.0,
        0.0,
        0.0,
    );
    println!("Earth position (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", earth_position.x().value());
    println!("  Y = {:.6} AU", earth_position.y().value());
    println!("  Z = {:.6} AU", earth_position.z().value());
    println!("  Distance from Sun = {:.6} AU\n", earth_position.distance().value());

    // Create a geocentric equatorial position (Moon at ~384,400 km)
    let moon_position = cartesian::position::Equatorial::<Kilometer>::new(
        300_000.0,
        200_000.0,
        100_000.0,
    );
    println!("Moon position (Geocentric Equatorial):");
    println!("  X = {:.1} km", moon_position.x().value());
    println!("  Y = {:.1} km", moon_position.y().value());
    println!("  Z = {:.1} km", moon_position.z().value());
    println!("  Distance from Earth = {:.1} km\n", moon_position.distance().value());

    // =========================================================================
    // 2. Spherical Coordinates
    // =========================================================================
    println!("2. SPHERICAL COORDINATES");
    println!("------------------------");

    // Create a star direction (Polaris approximately)
    let polaris = spherical::direction::Equatorial::new(
        Degrees::new(37.95),    // Right Ascension (converted to degrees)
        Degrees::new(89.26),    // Declination
    );
    println!("Polaris (Geocentric Equatorial Direction):");
    println!("  Right Ascension = {:.2}°", polaris.azimuth.value());
    println!("  Declination = {:.2}°\n", polaris.polar.value());

    // Create a position with distance (Betelgeuse at ~500 light-years)
    let betelgeuse_distance = 500.0 * 9.461e15 / 1.496e11; // Convert ly to AU
    let betelgeuse = spherical::position::ICRS::<AstronomicalUnit>::new(
        Degrees::new(88.79),    // RA
        Degrees::new(7.41),     // Dec
        betelgeuse_distance,
    );
    println!("Betelgeuse (Barycentric ICRS Position):");
    println!("  Right Ascension = {:.2}°", betelgeuse.azimuth.value());
    println!("  Declination = {:.2}°", betelgeuse.polar.value());
    println!("  Distance = {:.1} AU (~500 ly)\n", betelgeuse.distance.value());

    // =========================================================================
    // 3. Directions (Unit Vectors)
    // =========================================================================
    println!("3. DIRECTIONS (UNIT VECTORS)");
    println!("----------------------------");

    // Define an observer location for horizontal coordinates
    let observer = ObserverSite::new(
        Degrees::new(0.0),   // Longitude
        Degrees::new(51.5),  // Latitude
        0.0 * M,             // Height
    );

    // Directions are unitless (implicit radius = 1)
    let zenith = spherical::direction::Horizontal::with_site(
        observer,
        Degrees::new(90.0),     // Altitude (straight up)
        Degrees::new(0.0),      // Azimuth (North - doesn't matter for zenith)
    );
    println!("Zenith (Topocentric Horizontal from observer at 51.5°N):");
    println!("  Azimuth = {:.1}°", zenith.azimuth.value());
    println!("  Altitude = {:.1}°\n", zenith.polar.value());

    // Convert direction to position at a specific distance
    let cloud_distance = 5000.0 * KM;
    let cloud = zenith.position(cloud_distance);
    println!("Cloud at zenith, 5 km altitude:");
    println!("  Distance = {:.1} km\n", cloud.distance.value());

    // =========================================================================
    // 4. Cartesian <-> Spherical Conversion
    // =========================================================================
    println!("4. CARTESIAN ↔ SPHERICAL CONVERSION");
    println!("-----------------------------------");

    // Start with cartesian
    let cart_pos = cartesian::position::Equatorial::<AstronomicalUnit>::new(
        0.5,
        0.5,
        0.707,
    );
    println!("Cartesian position:");
    println!("  X = {:.3} AU", cart_pos.x().value());
    println!("  Y = {:.3} AU", cart_pos.y().value());
    println!("  Z = {:.3} AU", cart_pos.z().value());

    // Convert to spherical
    let sph_pos: spherical::position::Equatorial<AstronomicalUnit> = (&cart_pos).into();
    println!("\nConverted to Spherical:");
    println!("  RA = {:.2}°", sph_pos.azimuth.value());
    println!("  Dec = {:.2}°", sph_pos.polar.value());
    println!("  Distance = {:.3} AU", sph_pos.distance.value());

    // Convert back to cartesian
    let cart_pos_back: cartesian::position::Equatorial<AstronomicalUnit> = (&sph_pos).into();
    println!("\nConverted back to Cartesian:");
    println!("  X = {:.3} AU", cart_pos_back.x().value());
    println!("  Y = {:.3} AU", cart_pos_back.y().value());
    println!("  Z = {:.3} AU\n", cart_pos_back.z().value());

    // =========================================================================
    // 5. Type Safety
    // =========================================================================
    println!("5. TYPE SAFETY");
    println!("--------------");

    // Different coordinate types are incompatible
    let helio_pos = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    let geo_pos = cartesian::position::Equatorial::<AstronomicalUnit>::new(0.0, 1.0, 0.0);

    println!("Type-safe coordinates prevent mixing incompatible systems:");
    println!("  Heliocentric Ecliptic: {}", helio_pos);
    println!("  Geocentric Equatorial: {}", geo_pos);
    println!("\n  Cannot directly compute distance between them!");
    println!("  (Must transform to same center/frame first)\n");

    // But operations within the same type are allowed
    let pos1 = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    let pos2 = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.5, 0.0, 0.0);
    let distance = pos1.distance_to(&pos2);
    println!("Distance between two Heliocentric Ecliptic positions:");
    println!("  {:.3} AU\n", distance.value());

    // =========================================================================
    // 6. Different Centers and Frames
    // =========================================================================
    println!("6. CENTERS AND FRAMES");
    println!("---------------------");

    println!("Reference Centers:");
    println!("  Barycentric:  {}", centers::Barycentric::center_name());
    println!("  Heliocentric: {}", centers::Heliocentric::center_name());
    println!("  Geocentric:   {}", centers::Geocentric::center_name());
    println!("  Topocentric:  {}", centers::Topocentric::center_name());
    println!("  Bodycentric:  {}\n", centers::Bodycentric::center_name());

    println!("Reference Frames:");
    println!("  Ecliptic:   {}", frames::Ecliptic::frame_name());
    println!("  Equatorial: {}", frames::Equatorial::frame_name());
    println!("  Horizontal: {}", frames::Horizontal::frame_name());
    println!("  ICRS:       {}", frames::ICRS::frame_name());
    println!("  ECEF:       {}\n", frames::ECEF::frame_name());

    println!("=== Example Complete ===");
}
