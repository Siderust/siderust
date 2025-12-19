//! Observer-Based Coordinates Example
//!
//! This example demonstrates topocentric (observer-based) coordinates:
//! - Defining observer locations
//! - Horizontal coordinate system (altitude/azimuth)
//! - Converting between geocentric and topocentric

use qtty::*;
use siderust::astro::JulianDate;
use siderust::coordinates::cartesian::position::Equatorial;
use siderust::coordinates::centers::{Geocentric, ObserverSite};
use siderust::coordinates::spherical;

fn main() {
    println!("=== Observer-Based Coordinates Example ===\n");

    let jd = JulianDate::J2000;

    // =========================================================================
    // 1. Defining Observer Locations
    // =========================================================================
    println!("1. DEFINING OBSERVER LOCATIONS");
    println!("------------------------------");

    // Greenwich Observatory (Prime Meridian)
    let greenwich = ObserverSite::new(
        Degrees::new(0.0),     // Longitude: 0° (by definition)
        Degrees::new(51.4769), // Latitude: 51.4769° N
        0.0 * M,               // Height: sea level
    );

    println!("Greenwich Observatory:");
    println!("  Longitude: {:.4}° E", greenwich.lon);
    println!("  Latitude:  {:.4}° N", greenwich.lat);
    println!("  Height:    {:.0} m\n", greenwich.height);

    // Roque de los Muchachos Observatory (La Palma, Canary Islands)
    let la_palma = ObserverSite::new(
        Degrees::new(-17.8947), // Longitude: 17.8947° W
        Degrees::new(28.7606),  // Latitude: 28.7606° N
        2396.0 * M,             // Height: 2,396 m
    );

    println!("Roque de los Muchachos Observatory:");
    println!("  Longitude: {:.4}° W", -la_palma.lon);
    println!("  Latitude:  {:.4}° N", la_palma.lat);
    println!("  Height:    {:.0} m\n", la_palma.height);

    // Mauna Kea Observatory (Hawaii)
    let mauna_kea = ObserverSite::new(
        Degrees::new(-155.4783), // Longitude: 155.4783° W
        Degrees::new(19.8207),   // Latitude: 19.8207° N
        4207.0 * M,              // Height: 4,207 m
    );

    println!("Mauna Kea Observatory:");
    println!("  Longitude: {:.4}° W", -mauna_kea.lon);
    println!("  Latitude:  {:.4}° N", mauna_kea.lat);
    println!("  Height:    {:.0} m\n", mauna_kea.height);

    // Sydney Opera House (Australia)
    let sydney = ObserverSite::new(
        Degrees::new(151.2153), // Longitude: 151.2153° E
        Degrees::new(-33.8568), // Latitude: 33.8568° S (negative for south)
        0.0 * M,                // Height: sea level
    );

    println!("Sydney, Australia:");
    println!("  Longitude: {:.4}° E", sydney.lon);
    println!("  Latitude:  {:.4}° S", -sydney.lat);
    println!("  Height:    {:.0} m\n", sydney.height);

    // =========================================================================
    // 2. Horizontal Coordinate System
    // =========================================================================
    println!("2. HORIZONTAL COORDINATE SYSTEM");
    println!("-------------------------------");

    println!("The horizontal coordinate system is local to each observer:");
    println!("- Azimuth: Angle measured eastward from north (0° = North, 90° = East)");
    println!("- Altitude: Angle above the horizon (0° = horizon, 90° = zenith)\n");

    // Zenith (straight up) - a pure direction in horizontal frame
    // Note: Directions are now frame-only (no observer/center parameter)
    // The observer site is used for Positions, not for Directions
    let zenith = spherical::direction::Horizontal::new_horizontal(
        Degrees::new(90.0), // Altitude (straight up)
        Degrees::new(0.0),  // Azimuth (doesn't matter for zenith)
    );
    println!("Zenith direction:");
    println!("  Azimuth:  {:.1}°", zenith.az());
    println!("  Altitude: {:.1}° (straight up)\n", zenith.alt());

    // North horizon
    let north = spherical::direction::Horizontal::new_horizontal(
        Degrees::new(0.0), // On horizon
        Degrees::new(0.0), // North
    );
    println!("North point on horizon:");
    println!("  Azimuth:  {:.1}° (North)", north.az());
    println!("  Altitude: {:.1}° (horizon)\n", north.alt());

    // East, 30° above horizon
    let east_30 = spherical::direction::Horizontal::new_horizontal(
        Degrees::new(30.0), // 30° above horizon
        Degrees::new(90.0), // East
    );
    println!("30° above eastern horizon:");
    println!("  Azimuth:  {:.1}° (East)", east_30.az());
    println!("  Altitude: {:.1}°\n", east_30.alt());

    // =========================================================================
    // 3. Topocentric Positions (with distance)
    // =========================================================================
    println!("3. TOPOCENTRIC POSITIONS");
    println!("------------------------");

    // Create a geocentric position (e.g., a satellite)
    let satellite_geo: Equatorial<Kilometer, Geocentric> = Equatorial::new(
        5000.0, // X
        3000.0, // Y
        2000.0, // Z
    );

    println!("Satellite (Geocentric Equatorial):");
    println!("  X = {:.1} km", satellite_geo.x());
    println!("  Y = {:.1} km", satellite_geo.y());
    println!("  Z = {:.1} km", satellite_geo.z());
    println!(
        "  Distance from Earth center: {:.1} km\n",
        satellite_geo.distance()
    );

    // Transform to topocentric (Greenwich)
    let satellite_topo = satellite_geo.to_topocentric(greenwich, jd);

    println!("Satellite (Topocentric from Greenwich):");
    println!("  X = {:.1} km", satellite_topo.x());
    println!("  Y = {:.1} km", satellite_topo.y());
    println!("  Z = {:.1} km", satellite_topo.z());
    println!("  Note: Currently no parallax correction applied\n");

    // =========================================================================
    // 4. Multiple Observer Perspectives
    // =========================================================================
    println!("4. MULTIPLE OBSERVER PERSPECTIVES");
    println!("---------------------------------");

    // Same object viewed from different locations
    let object_geo: Equatorial<Kilometer, Geocentric> = Equatorial::new(10000.0, 0.0, 0.0);

    println!("Object position (Geocentric):");
    println!("  Distance: {:.0} km\n", object_geo.distance());

    println!("Same object from different observatories:");

    let from_greenwich = object_geo.to_topocentric(greenwich, jd);
    println!("  From Greenwich:");
    println!("    Observer lat: {:.2}° N", greenwich.lat);
    println!("    Distance: {:.0} km", from_greenwich.distance());

    let from_la_palma = object_geo.to_topocentric(la_palma, jd);
    println!("  From La Palma:");
    println!("    Observer lat: {:.2}° N", la_palma.lat);
    println!("    Distance: {:.0} km", from_la_palma.distance());

    let from_sydney = object_geo.to_topocentric(sydney, jd);
    println!("  From Sydney:");
    println!("    Observer lat: {:.2}° S", -sydney.lat);
    println!("    Distance: {:.0} km\n", from_sydney.distance());

    // =========================================================================
    // 5. Altitude-Dependent Phenomena
    // =========================================================================
    println!("5. ALTITUDE-DEPENDENT OBSERVATIONS");
    println!("----------------------------------");

    println!("Observer altitude affects:");
    println!("- Atmospheric pressure and refraction");
    println!("- Distance to visible horizon");
    println!("- Time of sunrise/sunset\n");

    let sea_level = ObserverSite::new(Degrees::new(0.0), Degrees::new(0.0), 0.0 * M);
    println!("Sea level observer:");
    println!("  Height: {:.0} m", sea_level.height);

    let mountain_top = ObserverSite::new(Degrees::new(0.0), Degrees::new(0.0), 4000.0 * M);
    println!("Mountain top observer:");
    println!("  Height: {:.0} m", mountain_top.height);

    // Approximate horizon distance formula: d ≈ 3.57 × √h (km)
    let horizon_sea = 3.57 * (sea_level.height.value() / 1000.0).sqrt();
    let horizon_mountain = 3.57 * (mountain_top.height.value() / 1000.0).sqrt();

    println!("\nApproximate distance to horizon:");
    println!("  From sea level: {:.1} km", horizon_sea);
    println!("  From 4000m: {:.1} km", horizon_mountain);
    println!("  Difference: {:.1} km\n", horizon_mountain - horizon_sea);

    // =========================================================================
    // 6. Converting Between Coordinate Systems
    // =========================================================================
    println!("6. COORDINATE SYSTEM CONVERSIONS");
    println!("--------------------------------");

    // Geocentric equatorial position
    let pos_geo = Equatorial::<Kilometer, Geocentric>::new(7000.0, 0.0, 0.0);
    println!("Position (Geocentric Equatorial):");
    println!("  Distance: {:.0} km", pos_geo.distance());

    // Convert to topocentric
    let pos_topo = pos_geo.to_topocentric(greenwich, jd);
    println!("\nPosition (Topocentric Equatorial from Greenwich):");
    println!("  Distance: {:.0} km", pos_topo.distance());

    // Convert topocentric to horizontal would require additional transforms
    // (currently requires sidereal time calculation for full implementation)

    println!("\n=== Example Complete ===");
    println!("\nKey Takeaways:");
    println!("- ObserverSite defines observer location (lon, lat, height)");
    println!("- Topocentric coordinates are relative to observer's location");
    println!("- Horizontal coordinates (altitude/azimuth) are observer-specific");
    println!("- Observer altitude affects horizon distance and observations");
    println!("- Different observers see different topocentric coordinates");
}
