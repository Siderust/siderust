// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orbit model overview.
//!
//! Shows:
//! - `KeplerianOrbit` as the canonical elliptic element set
//! - `MeanMotionOrbit` for mean-motion-driven propagation
//! - `ConicOrbit` for elliptic + hyperbolic classification/propagation
//! - `PreparedOrbit` for repeated fast propagation
//!
//! Run with: `cargo run --example 15_orbit_models`

use qtty::*;
use siderust::time::JulianDate;
use siderust::{ConicOrbit, KeplerianOrbit, MeanMotionOrbit, PreparedOrbit};
use std::convert::TryFrom;

fn main() {
    let jd = JulianDate::new(2_458_850.0);

    println!("=== Orbit Models ===\n");
    println!("Epoch: {jd}\n");

    let kepler = KeplerianOrbit::new(
        1.0 * AU,
        0.0167,
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(102.9),
        Degrees::new(100.0),
        JulianDate::J2000,
    );
    let kepler_pos = kepler.kepler_position(jd);
    println!("1. KeplerianOrbit");
    println!(
        "   heliocentric position = ({:.12}, {:.12}, {:.12}) AU",
        kepler_pos.x().value(),
        kepler_pos.y().value(),
        kepler_pos.z().value()
    );
    println!("   radius = {:.12} AU\n", kepler_pos.distance().value());

    let mean_motion = MeanMotionOrbit::try_new(
        1.0 * AU,
        0.0167,
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(102.9),
        0.9856,
        JulianDate::J2000,
    )
    .expect("valid mean-motion orbit");
    let mm_pos = mean_motion
        .position_at(jd)
        .expect("mean-motion propagation should succeed");
    println!("2. MeanMotionOrbit");
    println!(
        "   heliocentric position = ({:.12}, {:.12}, {:.12}) AU",
        mm_pos.x().value(),
        mm_pos.y().value(),
        mm_pos.z().value()
    );
    println!("   radius = {:.12} AU\n", mm_pos.distance().value());

    let halley_like = ConicOrbit::try_new(
        0.586 * AU,
        0.967,
        Degrees::new(162.3),
        Degrees::new(58.4),
        Degrees::new(111.3),
        Degrees::new(38.4),
        JulianDate::J2000,
    )
    .expect("valid elliptic conic");
    let hyperbolic = ConicOrbit::try_new(
        0.255 * AU,
        1.2,
        Degrees::new(122.7),
        Degrees::new(24.6),
        Degrees::new(241.8),
        Degrees::new(12.0),
        JulianDate::J2000,
    )
    .expect("valid hyperbolic conic");
    let elliptic_pos = halley_like
        .position_at(jd)
        .expect("elliptic conic propagation should succeed");
    let hyperbolic_pos = hyperbolic
        .position_at(jd)
        .expect("hyperbolic conic propagation should succeed");
    println!("3. ConicOrbit");
    println!("   halley-like kind = {:?}", halley_like.kind());
    println!("   hyperbolic kind  = {:?}", hyperbolic.kind());
    println!(
        "   elliptic radius  = {:.12} AU | hyperbolic radius = {:.12} AU\n",
        elliptic_pos.distance().value(),
        hyperbolic_pos.distance().value()
    );

    let prepared = PreparedOrbit::try_from(kepler).expect("elliptic orbit should prepare");
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
