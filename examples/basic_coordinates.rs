//! Basic Coordinates Example
//!
//! This example demonstrates the fundamental concepts of the siderust coordinate system:
//! - Creating positions and directions
//! - Different reference frames and centers
//! - Cartesian vs Spherical representations
//! - Type safety and units

use qtty::*;
use siderust::coordinates::cartesian;
use siderust::coordinates::centers::{self, ReferenceCenter};
use siderust::coordinates::frames::{self, ReferenceFrame};
use siderust::coordinates::spherical;

fn main() {
    println!("=== Siderust Basic Coordinates Example ===\n");

    // =========================================================================
    // 1. Cartesian Coordinates
    // =========================================================================
    println!("1. CARTESIAN COORDINATES");
    println!("------------------------");

    // Create a heliocentric ecliptic position (1 AU along X-axis)
    let earth_position = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    println!("Earth position (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", earth_position.x());
    println!("  Y = {:.6} AU", earth_position.y());
    println!("  Z = {:.6} AU", earth_position.z());
    println!(
        "  Distance from Sun = {:.6} AU\n",
        earth_position.distance()
    );

    // Create a geocentric equatorial position (Moon at ~384,400 km)
    let moon_position = cartesian::position::EquatorialMeanJ2000::<Kilometer>::new(
        300_000.0, 200_000.0, 100_000.0,
    );
    println!("Moon position (Geocentric EquatorialMeanJ2000):");
    println!("  X = {:.1} km", moon_position.x());
    println!("  Y = {:.1} km", moon_position.y());
    println!("  Z = {:.1} km", moon_position.z());
    println!(
        "  Distance from Earth = {:.1} km\n",
        moon_position.distance()
    );

    // =========================================================================
    // 2. Spherical Coordinates
    // =========================================================================
    println!("2. SPHERICAL COORDINATES");
    println!("------------------------");

    // Create a star direction (Polaris approximately)
    let polaris = spherical::direction::EquatorialMeanJ2000::new(
        Degrees::new(37.95), // Right Ascension (converted to degrees)
        Degrees::new(89.26), // Declination
    );
    println!("Polaris (Geocentric EquatorialMeanJ2000 Direction):");
    println!("  Right Ascension = {:.2}°", polaris.azimuth);
    println!("  Declination = {:.2}°\n", polaris.polar);

    // Create a position with distance (Betelgeuse at ~500 light-years)
    let betelgeuse_distance = 500.0 * 9.461e15 / 1.496e11; // Convert ly to AU
    let betelgeuse = spherical::position::ICRS::<AstronomicalUnit>::new(
        Degrees::new(88.79), // RA
        Degrees::new(7.41),  // Dec
        betelgeuse_distance,
    );
    println!("Betelgeuse (Barycentric ICRS Position):");
    println!("  Right Ascension = {:.2}°", betelgeuse.azimuth);
    println!("  Declination = {:.2}°", betelgeuse.polar);
    println!("  Distance = {:.1} AU (~500 ly)\n", betelgeuse.distance);

    // =========================================================================
    // 3. Directions (Unit Vectors)
    // =========================================================================
    println!("3. DIRECTIONS (UNIT VECTORS)");
    println!("----------------------------");

    // Directions are unitless (implicit radius = 1) and frame-only (no center)
    // Note: Directions don't carry observer site - they're pure directions
    let zenith = spherical::direction::Horizontal::new(
        Degrees::new(90.0), // Altitude (straight up)
        Degrees::new(0.0),  // Azimuth (North - doesn't matter for zenith)
    );
    println!("Zenith direction (Horizontal frame):");
    println!("  Altitude = {:.1}°", zenith.alt());
    println!("  Azimuth = {:.1}°\n", zenith.az());

    // Convert direction to position at a specific distance
    // Using Geocentric since it has simple Params = ()
    let cloud_distance = 5000.0 * KM;
    use siderust::coordinates::centers::Geocentric;
    let cloud = zenith.position::<Geocentric, _>(cloud_distance);
    println!("Cloud at zenith, 5 km altitude (relative to geocenter):");
    println!("  Distance = {:.1} km\n", cloud.distance);

    // =========================================================================
    // 4. Cartesian <-> Spherical Conversion
    // =========================================================================
    println!("4. CARTESIAN ↔ SPHERICAL CONVERSION");
    println!("-----------------------------------");

    // Start with cartesian
    let cart_pos =
        cartesian::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(0.5, 0.5, 0.707);
    println!("Cartesian position:");
    println!("  X = {:.3} AU", cart_pos.x());
    println!("  Y = {:.3} AU", cart_pos.y());
    println!("  Z = {:.3} AU", cart_pos.z());

    // Convert to spherical
    let sph_pos =
        spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::from_cartesian(&cart_pos);
    println!("\nConverted to Spherical:");
    println!("  RA = {:.2}°", sph_pos.azimuth);
    println!("  Dec = {:.2}°", sph_pos.polar);
    println!("  Distance = {:.3} AU", sph_pos.distance);

    // Convert back to cartesian
    let cart_pos_back =
        cartesian::position::EquatorialMeanJ2000::<AstronomicalUnit>::from_spherical(&sph_pos);
    println!("\nConverted back to Cartesian:");
    println!("  X = {:.3} AU", cart_pos_back.x());
    println!("  Y = {:.3} AU", cart_pos_back.y());
    println!("  Z = {:.3} AU\n", cart_pos_back.z());

    // =========================================================================
    // 5. Type Safety
    // =========================================================================
    println!("5. TYPE SAFETY");
    println!("--------------");

    // Different coordinate types are incompatible
    let helio_pos = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    let geo_pos =
        cartesian::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(0.0, 1.0, 0.0);

    println!("Type-safe coordinates prevent mixing incompatible systems:");
    println!("  Heliocentric Ecliptic: {}", helio_pos);
    println!("  Geocentric EquatorialMeanJ2000: {}", geo_pos);
    println!("\n  Cannot directly compute distance between them!");
    println!("  (Must transform to same center/frame first)\n");

    // But operations within the same type are allowed
    let pos1 = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    let pos2 = cartesian::position::Ecliptic::<AstronomicalUnit>::new(1.5, 0.0, 0.0);
    let distance = pos1.distance_to(&pos2);
    println!("Distance between two Heliocentric Ecliptic positions:");
    println!("  {:.3} AU\n", distance);

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
    println!(
        "  EquatorialMeanJ2000: {}",
        frames::EquatorialMeanJ2000::frame_name()
    );
    println!("  Horizontal: {}", frames::Horizontal::frame_name());
    println!("  ICRS:       {}", frames::ICRS::frame_name());
    println!("  ECEF:       {}\n", frames::ECEF::frame_name());

    println!("=== Example Complete ===");
}
