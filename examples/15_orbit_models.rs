// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orbit Models Example
//!
//! Demonstrates the four orbit propagation models: Keplerian, MeanMotion,
//! Conic, and Prepared.
//!
//! Run with: `cargo run --example 15_orbit_models`

use siderust::qtty::angular_rate::AngularRate;
use siderust::qtty::unit::{Day, Degree};
use siderust::qtty::AstronomicalUnits;
use siderust::time::JulianDate;
use siderust::{ConicKind, ConicOrbit, KeplerianOrbit, MeanMotionOrbit, PreparedOrbit};

fn main() {
    let jd = JulianDate::new(2458850.0);

    println!("=== Orbit Models ===\n");
    println!("Epoch: {:.1}\n", jd.jd_value());

    // 1. KeplerianOrbit
    let kepler = KeplerianOrbit::new(
        AstronomicalUnits::new(1.0),
        0.0167,
        siderust::qtty::Degrees::new(0.0),
        siderust::qtty::Degrees::new(0.0),
        siderust::qtty::Degrees::new(102.9),
        siderust::qtty::Degrees::new(100.0),
        JulianDate::new(2451545.0),
    );
    let kepler_pos = kepler.kepler_position(jd);
    println!("1. KeplerianOrbit");
    println!(
        "   heliocentric position = ({:.12}, {:.12}, {:.12}) AU",
        kepler_pos.x().value(),
        kepler_pos.y().value(),
        kepler_pos.z().value(),
    );
    println!("   radius = {:.12} AU\n", kepler_pos.distance().value());

    // 2. MeanMotionOrbit
    let mean_motion = MeanMotionOrbit::try_new(
        AstronomicalUnits::new(1.0),
        0.0167,
        siderust::qtty::Degrees::new(0.0),
        siderust::qtty::Degrees::new(0.0),
        siderust::qtty::Degrees::new(102.9),
        AngularRate::<Degree, Day>::new(0.9856),
        JulianDate::new(2451545.0),
    )
    .expect("valid mean-motion orbit");
    let mean_motion_pos = mean_motion.position_at(jd).expect("valid propagation");
    println!("2. MeanMotionOrbit");
    println!(
        "   heliocentric position = ({:.12}, {:.12}, {:.12}) AU",
        mean_motion_pos.x().value(),
        mean_motion_pos.y().value(),
        mean_motion_pos.z().value(),
    );
    println!(
        "   radius = {:.12} AU\n",
        mean_motion_pos.distance().value()
    );

    // 3. ConicOrbit
    let halley_like = ConicOrbit::try_new(
        AstronomicalUnits::new(0.586),
        0.967,
        siderust::qtty::Degrees::new(162.3),
        siderust::qtty::Degrees::new(58.4),
        siderust::qtty::Degrees::new(111.3),
        siderust::qtty::Degrees::new(38.4),
        JulianDate::new(2451545.0),
    )
    .expect("valid halley-like conic orbit");
    let hyperbolic = ConicOrbit::try_new(
        AstronomicalUnits::new(0.255),
        1.2,
        siderust::qtty::Degrees::new(122.7),
        siderust::qtty::Degrees::new(24.6),
        siderust::qtty::Degrees::new(241.8),
        siderust::qtty::Degrees::new(12.0),
        JulianDate::new(2451545.0),
    )
    .expect("valid hyperbolic conic orbit");
    let halley_pos = halley_like
        .position_at(jd)
        .expect("valid halley propagation");
    let hyperbolic_pos = hyperbolic
        .position_at(jd)
        .expect("valid hyperbolic propagation");
    println!("3. ConicOrbit");
    println!(
        "   halley-like kind = {}",
        if halley_like.kind() == ConicKind::Elliptic {
            "Elliptic"
        } else {
            "Hyperbolic"
        }
    );
    println!(
        "   hyperbolic kind  = {}",
        if hyperbolic.kind() == ConicKind::Elliptic {
            "Elliptic"
        } else {
            "Hyperbolic"
        }
    );
    println!(
        "   elliptic radius  = {:.12} AU | hyperbolic radius = {:.12} AU\n",
        halley_pos.distance().value(),
        hyperbolic_pos.distance().value(),
    );

    // 4. PreparedOrbit
    let prepared = PreparedOrbit::try_from(kepler).expect("valid prepared orbit");
    let prepared_pos = prepared.position_at(jd);
    let delta = ((kepler_pos.x().value() - prepared_pos.x().value()).powi(2)
        + (kepler_pos.y().value() - prepared_pos.y().value()).powi(2)
        + (kepler_pos.z().value() - prepared_pos.z().value()).powi(2))
    .sqrt();
    println!("4. PreparedOrbit");
    println!(
        "   prepared radius = {:.12} AU",
        prepared_pos.distance().value()
    );
    println!("   chord delta vs direct Keplerian = {:.3e}", delta);
}
