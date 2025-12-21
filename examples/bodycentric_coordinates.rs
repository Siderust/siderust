//! Body-Centric Coordinates Example
//!
//! This example demonstrates the new body-centric coordinate system that allows
//! viewing positions from arbitrary orbiting bodies (satellites, planets, moons, etc.).

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::{Earth, Mars, Venus};
use siderust::coordinates::cartesian::position::{Ecliptic, Position};
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::centers::{Bodycentric, BodycentricParams, Geocentric, Heliocentric};
use siderust::coordinates::frames;
use siderust::coordinates::transform::centers::{FromBodycentricExt, ToBodycentricExt};
use siderust::coordinates::transform::TransformCenter;

fn main() {
    println!("=== Body-Centric Coordinates Example ===\n");

    let jd = JulianDate::J2000;
    println!("Reference time: J2000.0 (JD {})\n", jd.value());

    // =========================================================================
    // 1. Satellite-Centric Coordinates (ISS example)
    // =========================================================================
    println!("1. SATELLITE-CENTRIC COORDINATES");
    println!("--------------------------------");

    // Define an ISS-like orbit (low Earth orbit)
    let iss_orbit = Orbit::new(
        0.0000426 * AU,     // ~6,378 km (Earth radius) in AU
        0.001,              // Nearly circular
        Degrees::new(51.6), // ISS inclination
        Degrees::new(0.0),  // Longitude of ascending node
        Degrees::new(0.0),  // Argument of perihelion
        Degrees::new(0.0),  // Mean anomaly at epoch
        jd,                 // Epoch
    );

    println!("ISS Orbit:");
    println!(
        "  Semi-major axis: {:.6} AU (~{:.0} km)",
        iss_orbit.semi_major_axis,
        iss_orbit.semi_major_axis.value() * 149597870.7
    );
    println!("  Eccentricity: {:.4}", iss_orbit.eccentricity);
    println!("  Inclination: {:.2}°\n", iss_orbit.inclination);

    // Create parameters for ISS-centric coordinates
    let iss_params = BodycentricParams::geocentric(iss_orbit);

    // Get ISS position in geocentric coordinates
    let iss_pos_ecl = iss_orbit.kepler_position(jd);
    println!("ISS position (Geocentric Ecliptic):");
    println!("  X = {:.8} AU", iss_pos_ecl.x());
    println!("  Y = {:.8} AU", iss_pos_ecl.y());
    println!("  Z = {:.8} AU", iss_pos_ecl.z());
    println!(
        "  Distance from Earth: {:.8} AU (~{:.0} km)\n",
        iss_pos_ecl.distance(),
        iss_pos_ecl.distance().value() * 149597870.7
    );

    // Moon's approximate position (geocentric)
    let moon_geo: Position<Geocentric, frames::Ecliptic, Au> = Position::new(0.00257, 0.0, 0.0); // ~384,400 km

    println!("Moon position (Geocentric):");
    println!(
        "  Distance from Earth: {:.5} AU (~{:.0} km)\n",
        moon_geo.distance(),
        moon_geo.distance().value() * 149597870.7
    );

    // Transform to ISS-centric coordinates
    let moon_from_iss: Position<Bodycentric, frames::Ecliptic, Au> =
        moon_geo.to_bodycentric(iss_params, jd);

    println!("Moon as seen from ISS:");
    println!("  X = {:.6} AU", moon_from_iss.x());
    println!("  Y = {:.6} AU", moon_from_iss.y());
    println!("  Z = {:.6} AU", moon_from_iss.z());
    println!(
        "  Distance from ISS: {:.5} AU (~{:.0} km)\n",
        moon_from_iss.distance(),
        moon_from_iss.distance().value() * 149597870.7
    );

    // =========================================================================
    // 2. Planet-Centric Coordinates (Mars view of Earth)
    // =========================================================================
    println!("2. MARS-CENTRIC COORDINATES");
    println!("---------------------------");

    // Mars orbit (simplified - using real values would need more orbital elements)
    let mars_orbit = Orbit::new(
        1.524 * AU,          // Semi-major axis
        0.0934,              // Eccentricity
        Degrees::new(1.85),  // Inclination
        Degrees::new(49.56), // Longitude of ascending node
        Degrees::new(286.5), // Argument of perihelion
        Degrees::new(19.41), // Mean anomaly at J2000
        jd,
    );

    let mars_params = BodycentricParams::heliocentric(mars_orbit);

    println!("Mars Orbit:");
    println!("  Semi-major axis: {:.3} AU", mars_orbit.semi_major_axis);
    println!("  Eccentricity: {:.4}", mars_orbit.eccentricity);
    println!("  Inclination: {:.2}°\n", mars_orbit.inclination);

    // Get actual positions from VSOP87
    let earth_helio = *Earth::vsop87a(jd).get_position();
    let mars_helio = *Mars::vsop87a(jd).get_position();

    println!("Earth (Heliocentric):");
    println!("  Distance from Sun: {:.6} AU\n", earth_helio.distance());

    println!("Mars (Heliocentric):");
    println!("  Distance from Sun: {:.6} AU\n", mars_helio.distance());

    // View Earth from Mars (using approximate orbit)
    let earth_from_mars: Position<Bodycentric, frames::Ecliptic, Au> =
        earth_helio.to_bodycentric(mars_params, jd);

    println!("Earth as seen from Mars:");
    println!("  X = {:.6} AU", earth_from_mars.x());
    println!("  Y = {:.6} AU", earth_from_mars.y());
    println!("  Z = {:.6} AU", earth_from_mars.z());
    println!("  Distance: {:.6} AU\n", earth_from_mars.distance());

    // =========================================================================
    // 3. Venus-Centric View
    // =========================================================================
    println!("3. VENUS-CENTRIC COORDINATES");
    println!("----------------------------");

    let venus_orbit = Orbit::new(
        0.723 * AU,
        0.0067,
        Degrees::new(3.39),
        Degrees::new(76.68),
        Degrees::new(131.53),
        Degrees::new(50.42),
        jd,
    );

    let venus_params = BodycentricParams::heliocentric(venus_orbit);

    let venus_helio = *Venus::vsop87a(jd).get_position();
    println!("Venus (Heliocentric):");
    println!("  Distance from Sun: {:.6} AU\n", venus_helio.distance());

    // View Earth from Venus
    let earth_from_venus: Position<Bodycentric, frames::Ecliptic, Au> =
        earth_helio.to_bodycentric(venus_params, jd);

    println!("Earth as seen from Venus:");
    println!("  Distance: {:.6} AU\n", earth_from_venus.distance());

    // View Mars from Venus
    let mars_from_venus: Position<Bodycentric, frames::Ecliptic, Au> =
        mars_helio.to_bodycentric(venus_params, jd);

    println!("Mars as seen from Venus:");
    println!("  Distance: {:.6} AU\n", mars_from_venus.distance());

    // =========================================================================
    // 4. Round-Trip Transformation
    // =========================================================================
    println!("4. ROUND-TRIP TRANSFORMATION");
    println!("----------------------------");

    let original_earth = earth_helio;
    println!("Original Earth position (Heliocentric):");
    println!("  X = {:.10} AU", original_earth.x());
    println!("  Y = {:.10} AU", original_earth.y());
    println!("  Z = {:.10} AU\n", original_earth.z());

    // Transform to Mars-centric
    let earth_mars_centric: Position<Bodycentric, frames::Ecliptic, Au> =
        original_earth.to_bodycentric(mars_params, jd);

    println!("Transformed to Mars-centric:");
    println!(
        "  Distance from Mars: {:.6} AU\n",
        earth_mars_centric.distance()
    );

    // Transform back to geocentric (to match original center)
    let recovered_geo: Ecliptic<Au, Geocentric> = earth_mars_centric.to_geocentric(jd);

    // Then to heliocentric to compare
    let recovered_helio: Ecliptic<Au, Heliocentric> = recovered_geo.to_center(jd);

    println!("Recovered position (Heliocentric):");
    println!("  X = {:.10} AU", recovered_helio.x());
    println!("  Y = {:.10} AU", recovered_helio.y());
    println!("  Z = {:.10} AU\n", recovered_helio.z());

    let diff_x = (original_earth.x() - recovered_helio.x()).value();
    let diff_y = (original_earth.y() - recovered_helio.y()).value();
    let diff_z = (original_earth.z() - recovered_helio.z()).value();
    let diff = (diff_x * diff_x + diff_y * diff_y + diff_z * diff_z).sqrt();
    println!("Total difference: {:.2e} AU (should be small)\n", diff);

    // =========================================================================
    // 5. Directions in Different Reference Frames
    // =========================================================================
    println!("5. DIRECTIONS AS FREE VECTORS");
    println!("------------------------------");

    // Directions are now frame-only (no center parameter).
    // They cannot undergo center transformations - only frame rotations.
    let star_dir: Direction<frames::Equatorial> = Direction::normalize(0.707, 0.0, 0.707); // 45° declination, 0° RA

    println!("Star direction (Equatorial frame):");
    println!("  X = {:.3}", star_dir.x());
    println!("  Y = {:.3}", star_dir.y());
    println!("  Z = {:.3}\n", star_dir.z());

    // Directions are invariant under center translations (they're free vectors).
    // The direction to a distant star looks the same from Earth or ISS
    // because directions don't carry center information.
    println!("Note: Directions are free vectors - they represent 'which way'");
    println!("without reference to any origin. For distant stars, the direction");
    println!("appears the same from any location in the solar system.\n");

    println!("=== Example Complete ===");
    println!("\nKey Takeaways:");
    println!("- Body-centric coordinates work for any orbiting body");
    println!("- Satellite-centric: use OrbitReferenceCenter::Geocentric");
    println!("- Planet-centric: use OrbitReferenceCenter::Heliocentric");
    println!("- Directions are free vectors (no center, only frame)");
    println!("- Round-trip transformations preserve positions (within precision)");
}
