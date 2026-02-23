// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Solar System Bodies Example
//!
//! This example demonstrates working with solar system bodies:
//! - Computing planetary positions using VSOP87
//! - Planet-to-planet views
//! - Distances and relative positions

use qtty::*;
use siderust::bodies::solar_system::*;
use siderust::coordinates::cartesian::position::EclipticMeanJ2000;
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::transform::TransformCenter;
use siderust::time::JulianDate;

fn main() {
    println!("=== Solar System Bodies Example ===\n");

    let j2000 = JulianDate::J2000;
    // Current time using chrono
    let current = JulianDate::from_utc(chrono::Utc::now());

    println!("Times:");
    println!("  J2000.0: JD {:.1}", j2000);
    println!("  Current: JD {:.1}\n", current);

    // =========================================================================
    // 1. Planetary Positions at J2000
    // =========================================================================
    println!("1. PLANETARY POSITIONS AT J2000");
    println!("-------------------------------");

    // Mercury
    let mercury = Mercury::vsop87a(j2000);
    println!("Mercury (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", mercury.x());
    println!("  Y = {:.6} AU", mercury.y());
    println!("  Z = {:.6} AU", mercury.z());
    println!("  Distance from Sun: {:.6} AU\n", mercury.distance());

    // Venus
    let venus = Venus::vsop87a(j2000);
    println!("Venus (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", venus.distance());

    // Earth
    let earth = Earth::vsop87a(j2000);
    println!("Earth (Heliocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", earth.x());
    println!("  Y = {:.6} AU", earth.y());
    println!("  Z = {:.6} AU", earth.z());
    println!(
        "  Distance from Sun: {:.6} AU (should be ~1.0)\n",
        earth.distance()
    );

    // Mars
    let mars = Mars::vsop87a(j2000);
    println!("Mars (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", mars.distance());

    // Jupiter
    let jupiter = Jupiter::vsop87a(j2000);
    println!("Jupiter (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", jupiter.distance());

    // Saturn
    let saturn = Saturn::vsop87a(j2000);
    println!("Saturn (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", saturn.distance());

    // Uranus
    let uranus = Uranus::vsop87a(j2000);
    println!("Uranus (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", uranus.distance());

    // Neptune
    let neptune = Neptune::vsop87a(j2000);
    println!("Neptune (Heliocentric EclipticMeanJ2000):");
    println!("  Distance from Sun: {:.6} AU\n", neptune.distance());

    // =========================================================================
    // 2. Planet-to-Planet Distances
    // =========================================================================
    println!("2. INTER-PLANETARY DISTANCES AT J2000");
    println!("-------------------------------------");

    // Earth to Mars
    let earth_pos = earth;
    let mars_pos = mars;
    let earth_mars_distance = earth_pos.distance_to(&mars_pos);
    let earth_mars_distance_km = earth_mars_distance.to::<Kilometer>();
    println!("Earth to Mars:");
    println!(
        "  Distance: {} ({})\n",
        earth_mars_distance, earth_mars_distance_km
    );

    // Earth to Jupiter
    let jupiter_pos = jupiter;
    let earth_jupiter_distance = earth_pos.distance_to(&jupiter_pos);
    let earth_jupiter_distance_km = earth_jupiter_distance.to::<Kilometer>();
    println!("Earth to Jupiter:");
    println!(
        "  Distance: {} ({})\n",
        earth_jupiter_distance, earth_jupiter_distance_km
    );

    // Venus to Mercury (inner planets)
    let venus_pos = venus;
    let mercury_pos = mercury;
    let venus_mercury_distance = venus_pos.distance_to(&mercury_pos);
    let venus_mercury_distance_km = venus_mercury_distance.to::<Kilometer>();
    println!("Venus to Mercury:");
    println!(
        "  Distance: {} ({})\n",
        venus_mercury_distance, venus_mercury_distance_km
    );

    // =========================================================================
    // 3. Geocentric Positions
    // =========================================================================
    println!("3. GEOCENTRIC POSITIONS (AS SEEN FROM EARTH)");
    println!("--------------------------------------------");

    // Sun as seen from Earth (negative Earth position)
    let sun_geo: EclipticMeanJ2000<Au, Geocentric> = earth_pos.to_center(j2000);
    println!("Sun (Geocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", sun_geo.x());
    println!("  Y = {:.6} AU", sun_geo.y());
    println!("  Z = {:.6} AU", sun_geo.z());
    println!("  Distance from Earth: {:.6} AU\n", sun_geo.distance());

    // Mars as seen from Earth
    let mars_geo: EclipticMeanJ2000<Au, Geocentric> = mars_pos.to_center(j2000);
    println!("Mars (Geocentric EclipticMeanJ2000):");
    println!("  X = {:.6} AU", mars_geo.x());
    println!("  Y = {:.6} AU", mars_geo.y());
    println!("  Z = {:.6} AU", mars_geo.z());
    println!("  Distance from Earth: {:.6} AU\n", mars_geo.distance());

    // Venus as seen from Earth
    let venus_geo: EclipticMeanJ2000<Au, Geocentric> = venus_pos.to_center(j2000);
    println!("Venus (Geocentric EclipticMeanJ2000):");
    println!("  Distance from Earth: {:.6} AU\n", venus_geo.distance());

    // Jupiter as seen from Earth
    let jupiter_geo: EclipticMeanJ2000<Au, Geocentric> = jupiter_pos.to_center(j2000);
    println!("Jupiter (Geocentric EclipticMeanJ2000):");
    println!("  Distance from Earth: {:.6} AU\n", jupiter_geo.distance());

    // =========================================================================
    // 4. Time Evolution
    // =========================================================================
    println!("4. PLANETARY POSITIONS OVER TIME");
    println!("--------------------------------");

    println!("Mars position evolution:");

    let dates = [
        ("J2000.0", j2000),
        ("6 months later", JulianDate::new(j2000.value() + 182.5)),
        ("1 year later", JulianDate::new(j2000.value() + 365.25)),
    ];

    for (label, date) in dates.iter() {
        let mars_t = Mars::vsop87a(*date);
        println!("  {} (JD {:.1}):", label, date);
        println!("    Distance from Sun: {:.6} AU", mars_t.distance());

        let mars_geo_t: EclipticMeanJ2000<Au, Geocentric> = mars_t.to_center(*date);
        println!("    Distance from Earth: {:.6} AU", mars_geo_t.distance());
    }
    println!();

    // =========================================================================
    // 5. Barycentric vs Heliocentric
    // =========================================================================
    println!("5. BARYCENTRIC VS HELIOCENTRIC COORDINATES");
    println!("------------------------------------------");

    // Sun position in barycentric coordinates
    let sun_bary = Sun::vsop87e(j2000);
    println!("Sun (Barycentric EclipticMeanJ2000):");
    println!("  X = {:.8} AU", sun_bary.x());
    println!("  Y = {:.8} AU", sun_bary.y());
    println!("  Z = {:.8} AU", sun_bary.z());
    println!("  Distance from SSB: {:.8} AU", sun_bary.distance());
    println!("  (Sun is offset from solar system barycenter due to Jupiter)\n");

    // Earth in barycentric coordinates
    let earth_bary = Earth::vsop87e(j2000);
    println!("Earth (Barycentric EclipticMeanJ2000):");
    println!("  Distance from SSB: {:.6} AU\n", earth_bary.distance());

    // Jupiter effect on barycenter
    let jupiter_bary = Jupiter::vsop87e(j2000);
    println!("Jupiter (Barycentric EclipticMeanJ2000):");
    println!("  Distance from SSB: {:.6} AU", jupiter_bary.distance());
    println!("  (Jupiter's mass causes the barycenter to be outside the Sun)\n");

    // =========================================================================
    // 6. Orbital Periods
    // =========================================================================
    println!("6. ORBITAL INFORMATION");
    println!("----------------------");

    println!("Approximate orbital periods:");
    println!("  Mercury: ~88 days");
    println!("  Venus:   ~225 days");
    println!("  Earth:   ~365.25 days");
    println!("  Mars:    ~687 days");
    println!("  Jupiter: ~11.9 years");
    println!("  Saturn:  ~29.5 years");
    println!("  Uranus:  ~84 years");
    println!("  Neptune: ~165 years\n");

    // =========================================================================
    // 7. Current Positions
    // =========================================================================
    println!("7. CURRENT PLANETARY POSITIONS");
    println!("------------------------------");

    let now_earth = Earth::vsop87a(current);
    let now_mars = Mars::vsop87a(current);
    let now_jupiter = Jupiter::vsop87a(current);

    println!("Current positions (JD {:.1}):", current);
    println!("  Earth distance from Sun: {:.6} AU", now_earth.distance());
    println!("  Mars distance from Sun: {:.6} AU", now_mars.distance());
    println!(
        "  Jupiter distance from Sun: {:.6} AU\n",
        now_jupiter.distance()
    );

    let now_mars_geo: EclipticMeanJ2000<Au, Geocentric> = now_mars.to_center(current);
    let now_mars_geo_km = now_mars_geo.distance().to::<Kilometer>();
    println!(
        "Current Mars distance from Earth: {} ({})\n",
        now_mars_geo.distance(),
        now_mars_geo_km
    );

    println!("=== Example Complete ===");
    println!("\nKey Takeaways:");
    println!("- VSOP87 provides high-precision planetary positions");
    println!("- VSOP87A: heliocentric ecliptic coordinates");
    println!("- VSOP87E: barycentric ecliptic coordinates");
    println!("- Positions change over time - always specify Julian Date");
    println!("- Convert between heliocentric and geocentric with .to_center(jd)");
    println!("- Solar system barycenter is offset from Sun due to Jupiter");
}
