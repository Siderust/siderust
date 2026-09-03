// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Horizontal Coordinates Example
//!
//! Demonstrates equatorial ↔ horizontal transforms for spherical and cartesian
//! directions and positions, with round-trip checks.
//!
//! Positions use **mean-of-date + GMST** via [`.transform(jd)`]. Directions use
//! **true-of-date + GAST** via [`.to_horizontal`] / [`.to_equatorial`] (observer
//! site passed explicitly). See `doc/frames_and_centers.md` for frame details.
//!
//! Run with: `cargo run --example 19_horizontal_coordinates`
#![allow(clippy::print_stdout)]

use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Topocentric;
use siderust::coordinates::frames::{EquatorialMeanOfDate, Horizontal};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::{DirectionAstroExt, SphericalDirectionAstroExt, Transform};
use siderust::qtty::*;
use siderust::time::JulianDate;

fn main() {
    println!("=== Horizontal Coordinates Example ===\n");

    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let jd = JulianDate::new(2_459_015.5); // 2020-01-01 12:00 TT

    // Sirius apparent equatorial coordinates (degrees).
    let ra = Degrees::new(101.287);
    let dec = Degrees::new(-16.716);
    let distance = 384_400.0 * KM;

    println!("Observer: {}", ROQUE_DE_LOS_MUCHACHOS.name);
    println!("Epoch (TT): {jd:.6}\n");

    // =========================================================================
    // 1. Spherical positions (EquatorialMeanOfDate → Horizontal, GMST)
    // =========================================================================
    println!("1. SPHERICAL POSITIONS (GMST, .transform)");
    println!("-------------------------------------------");

    let eq_mod_dir = spherical::direction::EquatorialMeanOfDate::new(ra, dec);
    let eq_sph = eq_mod_dir.position_with_params::<Topocentric, Kilometer>(site, distance);

    println!("EquatorialMeanOfDate (topocentric):");
    println!(
        "  RA = {:.4}, Dec = {:.4}, dist = {}",
        eq_sph.ra(),
        eq_sph.dec(),
        eq_sph.distance
    );

    let hz_sph: spherical::Position<Topocentric, Horizontal, Kilometer> = eq_sph.transform(jd);
    let back_sph: spherical::Position<Topocentric, EquatorialMeanOfDate, Kilometer> =
        hz_sph.transform(jd);

    println!("Horizontal:");
    println!(
        "  Alt = {:.4}, Az = {:.4}, dist = {}",
        hz_sph.alt(),
        hz_sph.az(),
        hz_sph.distance
    );
    println!(
        "  Round-trip distance residual: {}\n",
        (eq_sph.distance - back_sph.distance).abs()
    );

    // =========================================================================
    // 2. Cartesian positions (EquatorialMeanOfDate → Horizontal, GMST)
    // =========================================================================
    println!("2. CARTESIAN POSITIONS (GMST, .transform)");
    println!("-------------------------------------------");

    let eq_cart = Position::<Topocentric, EquatorialMeanOfDate, Kilometer>::from_spherical(&eq_sph);

    let hz_cart: Position<Topocentric, Horizontal, Kilometer> = eq_cart.transform(jd);
    let back_cart: Position<Topocentric, EquatorialMeanOfDate, Kilometer> = hz_cart.transform(jd);

    println!(
        "EquatorialMeanOfDate (cartesian): dist = {}",
        eq_cart.distance()
    );
    println!(
        "Horizontal (cartesian): Alt = {:.4}, Az = {:.4}, dist = {}",
        hz_cart.to_spherical().alt(),
        hz_cart.to_spherical().az(),
        hz_cart.distance()
    );
    println!(
        "  Round-trip distance residual: {}\n",
        (eq_cart.distance() - back_cart.distance()).abs()
    );

    // =========================================================================
    // 3. Spherical directions (EquatorialTrueOfDate → Horizontal, GAST)
    // =========================================================================
    println!("3. SPHERICAL DIRECTIONS (GAST, .to_horizontal / .to_equatorial)");
    println!("----------------------------------------------------------------");

    let eq_tod = spherical::direction::EquatorialTrueOfDate::new(ra, dec);

    println!("EquatorialTrueOfDate:");
    println!("  RA = {:.4}, Dec = {:.4}", eq_tod.ra(), eq_tod.dec());

    let hz_tod = eq_tod.to_horizontal(&jd, &site);
    let back_tod = hz_tod.to_equatorial(&jd, &site);

    println!("Horizontal:");
    println!("  Alt = {:.4}, Az = {:.4}", hz_tod.alt(), hz_tod.az());
    println!(
        "  Round-trip angular separation: {}\n",
        eq_tod.angular_separation(&back_tod)
    );

    // =========================================================================
    // 4. Cartesian directions (EquatorialTrueOfDate → Horizontal, GAST)
    // =========================================================================
    println!("4. CARTESIAN DIRECTIONS (GAST, .to_horizontal / .to_equatorial)");
    println!("-----------------------------------------------------------------");

    let eq_cart_dir = eq_tod.to_cartesian();
    let hz_cart_dir = eq_cart_dir.to_horizontal(&jd, &site);
    let back_cart_dir = hz_cart_dir.to_equatorial(&jd, &site);

    let hz_cart_sph = hz_cart_dir.to_spherical();
    let back_cart_sph = back_cart_dir.to_spherical();

    println!("Horizontal (cartesian):");
    println!(
        "  Alt = {:.4}, Az = {:.4}",
        hz_cart_sph.alt(),
        hz_cart_sph.az()
    );
    println!(
        "  Round-trip angular separation: {}\n",
        eq_tod.angular_separation(&back_cart_sph)
    );

    println!("=== Example Complete ===");
}
