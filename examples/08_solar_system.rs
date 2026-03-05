// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Solar System + Planets Module Tour
//!
//! Run with: `cargo run --example 12_solar_system_example`

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::bodies::planets::{OrbitExt, Planet};
use siderust::bodies::solar_system::*;
use siderust::calculus::vsop87::VSOP87;
use siderust::coordinates::cartesian::position::EclipticMeanJ2000;
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::transform::TransformCenter;
use siderust::time::JulianDate;

fn main() {
    let jd = JulianDate::J2000;
    let now = JulianDate::from_utc(chrono::Utc::now());

    println!("=== Siderust Solar System Module Tour ===\n");
    println!("Epoch used for deterministic outputs: J2000 (JD {:.1})", jd);
    println!("Current epoch snapshot: JD {:.6}\n", now);

    section_catalog_overview();
    section_planet_constants_and_periods();
    section_vsop87_positions(jd);
    section_center_transforms(jd);
    section_moon_and_lagrange_points(jd);
    section_trait_dispatch(jd);
    section_planet_builder();
    section_current_snapshot(now);
}

fn section_catalog_overview() {
    println!("1) CATALOG OVERVIEW");
    println!("------------------");
    println!("Sun: {}", SOLAR_SYSTEM.sun.name);
    println!("Major planets: {}", SOLAR_SYSTEM.planets.len());
    println!("Dwarf planets: {}", SOLAR_SYSTEM.dwarf_planets.len());
    println!("Major moons: {}", SOLAR_SYSTEM.moons.len());
    println!("Lagrange points: {}\n", SOLAR_SYSTEM.lagrange_points.len());
}

fn section_planet_constants_and_periods() {
    println!("2) PLANET CONSTANTS + ORBITEXT::period()\n");

    let major_planets = [
        ("Mercury", &MERCURY),
        ("Venus", &VENUS),
        ("Earth", &EARTH),
        ("Mars", &MARS),
        ("Jupiter", &JUPITER),
        ("Saturn", &SATURN),
        ("Uranus", &URANUS),
        ("Neptune", &NEPTUNE),
    ];

    println!(
        "{:<8} {:>10} {:>10} {:>10}",
        "Planet", "a [AU]", "e", "Period"
    );
    println!("{}", "-".repeat(48));
    for (name, p) in major_planets {
        println!(
            "{:<8} {:>10.6} {:>10.6} {:>10.2}",
            name,
            p.orbit.semi_major_axis,
            p.orbit.eccentricity,
            p.orbit.period().to::<Day>()
        );
    }
    println!();
}

fn section_vsop87_positions(jd: JulianDate) {
    println!("3) VSOP87 EPHEMERIDES (HELIOCENTRIC + BARYCENTRIC)");
    println!("-----------------------------------------------");

    let earth_h = Earth::vsop87a(jd);
    let mars_h = Mars::vsop87a(jd);
    let earth_mars = earth_h.distance_to(&mars_h);

    println!("Earth heliocentric distance: {:.6}", earth_h.distance());
    println!("Mars heliocentric distance:  {:.6}", mars_h.distance());
    println!(
        "Earth-Mars separation:      {:.6} ({:.0}",
        earth_mars,
        earth_mars.to::<Kilometer>()
    );

    let sun_bary = Sun::vsop87e(jd);
    println!(
        "Sun barycentric offset from SSB: {:.8}\n",
        sun_bary.distance()
    );

    let (jupiter_pos, jupiter_vel) = Jupiter::vsop87e_pos_vel(jd);
    println!("Jupiter barycentric position+velocity at J2000:");
    println!("  r = {:.6}", jupiter_pos.x());
    println!("  v = {:.6}\n", jupiter_vel);
}

fn section_center_transforms(jd: JulianDate) {
    println!("4) CENTER TRANSFORMS (HELIOCENTRIC -> GEOCENTRIC)");
    println!("-----------------------------------------------");

    let mars_helio = Mars::vsop87a(jd);
    let mars_geo: EclipticMeanJ2000<AstronomicalUnit, Geocentric> = mars_helio.to_center(jd);

    println!(
        "Mars geocentric distance at J2000: {:.6}",
        mars_geo.distance()
    );
    println!(
        "Mars geocentric distance at J2000: {:.0}\n",
        mars_geo.distance().to::<Kilometer>()
    );
}

fn section_moon_and_lagrange_points(jd: JulianDate) {
    println!("5) MOON + LAGRANGE POINTS");
    println!("-------------------------");

    let moon_geo = Moon::get_geo_position::<Kilometer>(jd);
    println!(
        "Moon geocentric distance (ELP2000): {:.1} ({:.6})",
        moon_geo.distance(),
        moon_geo.distance().to::<AstronomicalUnit>()
    );

    println!("Lagrange points available in the catalog:");
    for lp in LAGRANGE_POINTS {
        println!(
            "  {:<12} in {:<10} -> lon={:>7.2} lat={:>6.2} r={:>5.2}",
            lp.name, lp.parent_system, lp.position.azimuth, lp.position.polar, lp.position.distance
        );
    }
    println!();
}

fn section_trait_dispatch(jd: JulianDate) {
    println!("6) TRAIT-BASED DISPATCH (calculus::vsop87::VSOP87)");
    println!("-------------------------------------------------");

    let dynamic_planets: [(&str, &dyn VSOP87); 4] = [
        ("Mercury", &Mercury),
        ("Venus", &Venus),
        ("Earth", &Earth),
        ("Mars", &Mars),
    ];

    for (name, planet) in dynamic_planets {
        let helio = planet.vsop87a(jd);
        let bary = planet.vsop87e(jd);
        println!(
            "{:<8} helio={:>8.5}  bary={:>8.5}",
            name,
            helio.distance(),
            bary.distance()
        );
    }
    println!();
}

fn section_planet_builder() {
    println!("7) planets::Planet BUILDER + OrbitExt");
    println!("-------------------------------------");

    let demo_world = Planet::builder()
        .mass(Kilograms::new(5.972e24 * 2.0))
        .radius(Kilometers::new(6371.0 * 1.3))
        .orbit(Orbit::new(
            AstronomicalUnits::new(1.4),
            0.07,
            Degrees::new(4.0),
            Degrees::new(120.0),
            Degrees::new(80.0),
            Degrees::new(10.0),
            JulianDate::J2000,
        ))
        .build();

    println!("Custom planet built at runtime:");
    println!("  mass   = {}", demo_world.mass);
    println!("  radius = {}", demo_world.radius);
    println!("  a      = {}", demo_world.orbit.semi_major_axis);
    println!(
        "  sidereal period = {:.2}\n",
        demo_world.orbit.period().to::<Day>()
    );
}

fn section_current_snapshot(now: JulianDate) {
    println!("8) CURRENT SNAPSHOT");
    println!("-------------------");

    let earth_now = Earth::vsop87a(now);
    let mars_now = Mars::vsop87a(now);
    let mars_geo_now: EclipticMeanJ2000<AstronomicalUnit, Geocentric> = mars_now.to_center(now);

    println!("Earth-Sun distance now: {:.6}", earth_now.distance());
    println!("Mars-Sun distance now:  {:.6}", mars_now.distance());
    println!(
        "Mars-Earth distance now: {:.6} ({:.0})",
        mars_geo_now.distance(),
        mars_geo_now.distance().to::<Kilometer>()
    );

    println!("\n=== End of example ===");
}
