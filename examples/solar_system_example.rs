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
use siderust::coordinates::cartesian::position::Ecliptic;
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
    println!("Mercury (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", mercury.get_position().x());
    println!("  Y = {:.6} AU", mercury.get_position().y());
    println!("  Z = {:.6} AU", mercury.get_position().z());
    println!(
        "  Distance from Sun: {:.6} AU\n",
        mercury.get_position().distance()
    );

    // Venus
    let venus = Venus::vsop87a(j2000);
    println!("Venus (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        venus.get_position().distance()
    );

    // Earth
    let earth = Earth::vsop87a(j2000);
    println!("Earth (Heliocentric Ecliptic):");
    println!("  X = {:.6} AU", earth.get_position().x());
    println!("  Y = {:.6} AU", earth.get_position().y());
    println!("  Z = {:.6} AU", earth.get_position().z());
    println!(
        "  Distance from Sun: {:.6} AU (should be ~1.0)\n",
        earth.get_position().distance()
    );

    // Mars
    let mars = Mars::vsop87a(j2000);
    println!("Mars (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        mars.get_position().distance()
    );

    // Jupiter
    let jupiter = Jupiter::vsop87a(j2000);
    println!("Jupiter (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        jupiter.get_position().distance()
    );

    // Saturn
    let saturn = Saturn::vsop87a(j2000);
    println!("Saturn (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        saturn.get_position().distance()
    );

    // Uranus
    let uranus = Uranus::vsop87a(j2000);
    println!("Uranus (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        uranus.get_position().distance()
    );

    // Neptune
    let neptune = Neptune::vsop87a(j2000);
    println!("Neptune (Heliocentric Ecliptic):");
    println!(
        "  Distance from Sun: {:.6} AU\n",
        neptune.get_position().distance()
    );

    // =========================================================================
    // 2. Planet-to-Planet Distances
    // =========================================================================
    println!("2. INTER-PLANETARY DISTANCES AT J2000");
    println!("-------------------------------------");

    // Earth to Mars
    let earth_pos = *earth.get_position();
    let mars_pos = *mars.get_position();
    let earth_mars_distance = earth_pos.distance_to(&mars_pos);
    let earth_mars_distance_km = earth_mars_distance.to::<Kilometer>();
    println!("Earth to Mars:");
    println!(
        "  Distance: {} ({})\n",
        earth_mars_distance,
        earth_mars_distance_km
    );

    // Earth to Jupiter
    let jupiter_pos = *jupiter.get_position();
    let earth_jupiter_distance = earth_pos.distance_to(&jupiter_pos);
    let earth_jupiter_distance_km = earth_jupiter_distance.to::<Kilometer>();
    println!("Earth to Jupiter:");
    println!(
        "  Distance: {} ({})\n",
        earth_jupiter_distance,
        earth_jupiter_distance_km
    );

    // Venus to Mercury (inner planets)
    let venus_pos = *venus.get_position();
    let mercury_pos = *mercury.get_position();
    let venus_mercury_distance = venus_pos.distance_to(&mercury_pos);
    let venus_mercury_distance_km = venus_mercury_distance.to::<Kilometer>();
    println!("Venus to Mercury:");
    println!(
        "  Distance: {} ({})\n",
        venus_mercury_distance,
        venus_mercury_distance_km
    );

    // =========================================================================
    // 3. Geocentric Positions
    // =========================================================================
    println!("3. GEOCENTRIC POSITIONS (AS SEEN FROM EARTH)");
    println!("--------------------------------------------");

    // Sun as seen from Earth (negative Earth position)
    let sun_geo: Ecliptic<Au, Geocentric> = earth_pos.to_center(j2000);
    println!("Sun (Geocentric Ecliptic):");
    println!("  X = {:.6} AU", sun_geo.x());
    println!("  Y = {:.6} AU", sun_geo.y());
    println!("  Z = {:.6} AU", sun_geo.z());
    println!("  Distance from Earth: {:.6} AU\n", sun_geo.distance());

    // Mars as seen from Earth
    let mars_geo: Ecliptic<Au, Geocentric> = mars_pos.to_center(j2000);
    println!("Mars (Geocentric Ecliptic):");
    println!("  X = {:.6} AU", mars_geo.x());
    println!("  Y = {:.6} AU", mars_geo.y());
    println!("  Z = {:.6} AU", mars_geo.z());
    println!("  Distance from Earth: {:.6} AU\n", mars_geo.distance());

    // Venus as seen from Earth
    let venus_geo: Ecliptic<Au, Geocentric> = venus_pos.to_center(j2000);
    println!("Venus (Geocentric Ecliptic):");
    println!("  Distance from Earth: {:.6} AU\n", venus_geo.distance());

    // Jupiter as seen from Earth
    let jupiter_geo: Ecliptic<Au, Geocentric> = jupiter_pos.to_center(j2000);
    println!("Jupiter (Geocentric Ecliptic):");
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
        println!(
            "    Distance from Sun: {:.6} AU",
            mars_t.get_position().distance()
        );

        let mars_geo_t: Ecliptic<Au, Geocentric> = mars_t.get_position().to_center(*date);
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
    println!("Sun (Barycentric Ecliptic):");
    println!("  X = {:.8} AU", sun_bary.get_position().x());
    println!("  Y = {:.8} AU", sun_bary.get_position().y());
    println!("  Z = {:.8} AU", sun_bary.get_position().z());
    println!(
        "  Distance from SSB: {:.8} AU",
        sun_bary.get_position().distance()
    );
    println!("  (Sun is offset from solar system barycenter due to Jupiter)\n");

    // Earth in barycentric coordinates
    let earth_bary = Earth::vsop87e(j2000);
    println!("Earth (Barycentric Ecliptic):");
    println!(
        "  Distance from SSB: {:.6} AU\n",
        earth_bary.get_position().distance()
    );

    // Jupiter effect on barycenter
    let jupiter_bary = Jupiter::vsop87e(j2000);
    println!("Jupiter (Barycentric Ecliptic):");
    println!(
        "  Distance from SSB: {:.6} AU",
        jupiter_bary.get_position().distance()
    );
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
    println!(
        "  Earth distance from Sun: {:.6} AU",
        now_earth.get_position().distance()
    );
    println!(
        "  Mars distance from Sun: {:.6} AU",
        now_mars.get_position().distance()
    );
    println!(
        "  Jupiter distance from Sun: {:.6} AU\n",
        now_jupiter.get_position().distance()
    );

    let now_mars_geo: Ecliptic<Au, Geocentric> = now_mars.get_position().to_center(current);
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
