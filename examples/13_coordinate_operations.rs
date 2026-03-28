// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Coordinate Operations Example
//!
//! Demonstrates what you can do with spherical/cartesian coordinate types:
//! angular separation (Vincenty formula), cross-validation against the
//! cartesian dot product, Euclidean 3D distance, and how type safety prevents
//! accidentally mixing incompatible coordinate systems.
//!
//! Run with: `cargo run --example 13_coordinate_operations`

use std::f64::consts::PI;

use qtty::*;
use siderust::coordinates::spherical;

fn main() {
    println!("=== Siderust Coordinate Operations Example ===\n");

    // =========================================================================
    // 1. Angular Separation, Spherical Directions
    // =========================================================================
    println!("1. ANGULAR SEPARATION (SPHERICAL DIRECTIONS)");
    println!("---------------------------------------------");

    // Two well-known stars using EquatorialMeanJ2000 directions.
    // Constructor: new(ra, dec), azimuth first, matching IAU naming.
    let polaris = spherical::direction::EquatorialMeanJ2000::new(
        Degrees::new(37.9546), // Right Ascension
        Degrees::new(89.2641), // Declination
    );
    let sirius = spherical::direction::EquatorialMeanJ2000::new(
        Degrees::new(101.2872), // Right Ascension
        Degrees::new(-16.7161), // Declination
    );

    let sep = polaris.angular_separation(&sirius);
    println!(
        "Polaris  (RA={:.4}°, Dec={:.4}°)",
        polaris.ra(),
        polaris.dec()
    );
    println!(
        "Sirius   (RA={:.4}°, Dec={:.4}°)",
        sirius.ra(),
        sirius.dec()
    );
    println!("  Angular separation = {:.4}°\n", sep);

    // Nearby pair, separation should be about 0.71°
    let star_a =
        spherical::direction::EquatorialMeanJ2000::new(Degrees::new(10.0), Degrees::new(30.0));
    let star_b =
        spherical::direction::EquatorialMeanJ2000::new(Degrees::new(10.5), Degrees::new(30.5));
    let close_sep = star_a.angular_separation(&star_b);
    println!("Nearby pair (0.5° apart in both RA and Dec):");
    println!("  Angular separation = {:.4}°\n", close_sep);

    // Self-separation must be exactly 0
    let self_sep = polaris.angular_separation(&polaris);
    println!(
        "Self-separation of Polaris = {:.6}°  (must be 0)\n",
        self_sep
    );

    // =========================================================================
    // 2. Cross-Validation: Spherical Vincenty vs. Cartesian Dot Product
    // =========================================================================
    println!("2. CROSS-VALIDATION: SPHERICAL vs. CARTESIAN");
    println!("---------------------------------------------");

    // Convert spherical directions to Cartesian unit vectors
    let polaris_cart = polaris.to_cartesian();
    let sirius_cart = sirius.to_cartesian();

    // angle_to returns radians
    let angle_rad = polaris_cart.angle_to(&sirius_cart);
    let angle_deg = angle_rad * 180.0 / PI;

    println!(
        "Polaris cartesian: ({:.4}, {:.4}, {:.4})",
        polaris_cart.x(),
        polaris_cart.y(),
        polaris_cart.z()
    );
    println!(
        "Sirius  cartesian: ({:.4}, {:.4}, {:.4})",
        sirius_cart.x(),
        sirius_cart.y(),
        sirius_cart.z()
    );
    println!(
        "  angle_to (Cartesian)         = {:.6} rad = {:.4}°",
        angle_rad, angle_deg
    );
    println!("  angular_separation (Vincenty) = {:.4}°", sep);
    println!(
        "  Difference                   = {:.2e}°  (near machine epsilon)\n",
        (angle_deg - sep.value()).abs()
    );

    // =========================================================================
    // 3. Dot Product and Perpendicularity
    // =========================================================================
    println!("3. DOT PRODUCT");
    println!("--------------");

    // North celestial pole and a point on the equator are perpendicular
    let north_pole =
        spherical::direction::EquatorialMeanJ2000::new(Degrees::new(0.0), Degrees::new(90.0))
            .to_cartesian();
    let equatorial =
        spherical::direction::EquatorialMeanJ2000::new(Degrees::new(0.0), Degrees::new(0.0))
            .to_cartesian();

    println!(
        "dot(North Pole, Equatorial point)   = {:.6}  (must be  0)",
        north_pole.dot(&equatorial)
    );

    // Anti-Polaris: opposite direction on the sky
    let anti_polaris = spherical::direction::EquatorialMeanJ2000::new(
        Degrees::new(polaris.ra().value() + 180.0),
        Degrees::new(-polaris.dec().value()),
    )
    .to_cartesian();
    println!(
        "dot(Polaris, anti-Polaris)          = {:.6}  (must be -1)\n",
        polaris_cart.dot(&anti_polaris)
    );

    // =========================================================================
    // 4. Euclidean Distance Between Spherical Positions
    // =========================================================================
    println!("4. EUCLIDEAN DISTANCE BETWEEN SPHERICAL POSITIONS");
    println!("--------------------------------------------------");

    // Approximate heliocentric ecliptic positions (lon, lat, distance)
    let earth = spherical::position::EclipticMeanJ2000::<AstronomicalUnit>::new(
        Degrees::new(100.0), // ecliptic longitude
        Degrees::new(0.0),   // ecliptic latitude
        1.0,                 // 1 AU
    );
    let mars = spherical::position::EclipticMeanJ2000::<AstronomicalUnit>::new(
        Degrees::new(200.0), // ecliptic longitude (~100° further)
        Degrees::new(2.0),   // ecliptic latitude
        1.524,               // ~1.524 AU
    );

    // Position is Copy, so we can pass by value freely.
    let eu_dist = earth.distance_to(&mars);
    let ang_sep_pos = earth.angular_separation(mars); // takes by value (Copy)

    println!("Earth (lon=100°, lat=0°, r=1.000 AU)");
    println!("Mars  (lon=200°, lat=2°, r=1.524 AU)");
    println!("  Euclidean 3D distance  = {:.4}", eu_dist);
    println!("  Angular separation     = {:.4}°\n", ang_sep_pos);

    // =========================================================================
    // 5. Same Operations from Cartesian Positions
    // =========================================================================
    println!("5. CARTESIAN POSITION OPERATIONS");
    println!("---------------------------------");

    let earth_cart = earth.to_cartesian();
    let mars_cart = mars.to_cartesian();

    let eu_dist_cart = earth_cart.distance_to(&mars_cart);
    println!("Earth cartesian: {earth_cart}");
    println!("Mars  cartesian: {mars_cart}");
    println!(
        "  Euclidean distance = {:.4}  (matches spherical)",
        eu_dist_cart
    );

    // Vector difference gives displacement vector
    let diff = mars_cart - earth_cart;
    println!("  Mars − Earth vector: {}", diff);
    println!("  |Mars − Earth|      = {:.4}\n", diff.magnitude());

    // =========================================================================
    // 6. Ecliptic Directions, Same API, Different Frame
    // =========================================================================
    println!("6. ECLIPTIC AND EQUATORIAL DIRECTIONS");
    println!("--------------------------------------");

    // EclipticMeanJ2000 uses lon/lat convention
    let vernal_equinox = spherical::direction::EclipticMeanJ2000::new(
        Degrees::new(0.0), // ecliptic longitude
        Degrees::new(0.0), // ecliptic latitude
    );
    let summer_solstice = spherical::direction::EclipticMeanJ2000::new(
        Degrees::new(90.0), // 90° along ecliptic
        Degrees::new(0.0),
    );

    let equinox_to_solstice = vernal_equinox.angular_separation(&summer_solstice);
    println!(
        "Vernal Equinox  → Summer Solstice angular sep = {:.4}°  (must be 90°)",
        equinox_to_solstice
    );

    let ecliptic_north = spherical::direction::EclipticMeanJ2000::new(
        Degrees::new(0.0),
        Degrees::new(90.0), // ecliptic north pole
    );
    let ecliptic_to_equinox = ecliptic_north.angular_separation(&vernal_equinox);
    println!(
        "Ecliptic North Pole → Vernal Equinox ang. sep = {:.4}°  (must be 90°)\n",
        ecliptic_to_equinox
    );

    println!("  Type safety note: EclipticMeanJ2000 and EquatorialMeanJ2000");
    println!("  directions are distinct types, the compiler prevents mixing them.");
    println!("  angular_separation only compiles for directions in the same frame.\n");

    println!("=== Example Complete ===");
}
