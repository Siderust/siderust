// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target examples.
//!
//! Shows how to use:
//! - `Trackable` for dynamic sky objects (Sun, planets, Moon, stars, ICRS directions)
//! - `Target<T>` / `CoordinateWithPM<T>` as timestamped coordinate snapshots
//! - optional proper motion for stellar targets
//! - target frame+center conversion through `From<&Target<_>>`
//!
//! Run with:
//! `cargo run --example 08_target`

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::astro::proper_motion::{set_proper_motion_since_j2000, ProperMotion};
use siderust::bodies::comet::HALLEY;
use siderust::bodies::solar_system::{Mars, Moon, Sun};
use siderust::bodies::{catalog, Satellite};
use siderust::coordinates::cartesian;
use siderust::coordinates::centers::{Geocentric, Heliocentric};
use siderust::coordinates::spherical::direction;
use siderust::targets::{Target, Trackable};
use siderust::time::JulianDate;

fn main() {
    let jd = JulianDate::J2000;
    let jd_next = jd + Days::new(1.0);

    println!("Target + Trackable examples");
    println!("===========================\n");

    section_trackable_objects(jd, jd_next);
    section_target_snapshots(jd, jd_next);
    section_target_with_proper_motion(jd);
    section_target_transform(jd);
}

fn section_trackable_objects(jd: JulianDate, jd_next: JulianDate) {
    println!("1) Trackable objects (ICRS, star, Sun, planet, Moon)");

    let fixed_icrs = direction::ICRS::new(Degrees::new(120.0), Degrees::new(22.5));
    let fixed_icrs_at_jd = fixed_icrs.track(jd);
    let fixed_icrs_at_next = fixed_icrs.track(jd_next);

    println!(
        "  ICRS direction is time-invariant: RA {:.3} -> {:.3}, Dec {:.3} -> {:.3}",
        fixed_icrs_at_jd.ra(),
        fixed_icrs_at_next.ra(),
        fixed_icrs_at_jd.dec(),
        fixed_icrs_at_next.dec()
    );

    let sirius_dir = catalog::SIRIUS.track(jd);
    println!(
        "  Sirius via Trackable: RA {:.3}, Dec {:.3}",
        sirius_dir.ra(),
        sirius_dir.dec()
    );

    let sun = Sun.track(jd);
    let mars = Mars.track(jd);
    let moon = Moon.track(jd);

    println!("  Sun barycentric distance: {:.6}", sun.position.distance());
    println!(
        "  Mars barycentric distance: {:.6}",
        mars.position.distance()
    );
    println!("  Moon geocentric distance: {:.1}\n", moon.distance());
}

fn section_target_snapshots(jd: JulianDate, jd_next: JulianDate) {
    println!("2) Target snapshots for arbitrary sky objects");

    // Planet snapshot target (heliocentric ecliptic cartesian)
    let mut mars_target = Target::new_static(Mars::vsop87a(jd), jd);
    println!(
        "  Mars target at JD {:.1}: r = {:.6}",
        mars_target.time,
        mars_target.position.distance()
    );

    // Update target with a new ephemeris sample at the next epoch.
    mars_target.update(Mars::vsop87a(jd_next), jd_next);
    println!(
        "  Mars target updated to JD {:.1}: r = {:.6}",
        mars_target.time,
        mars_target.position.distance()
    );

    // Comet snapshot target (orbit propagated with Kepler helper).
    let halley_target = Target::new_static(HALLEY.orbit.kepler_position(jd), jd);
    println!(
        "  Halley target at JD {:.1}: r = {:.6}",
        halley_target.time,
        halley_target.position.distance()
    );

    // Satellite-like custom object: propagate its Orbit and wrap in Target.
    let demo_satellite = Satellite::new(
        "DemoSat",
        Kilograms::new(1_200.0),
        Kilometers::new(1.4),
        Orbit::new(
            AstronomicalUnits::new(1.0002),
            0.001,
            Degrees::new(0.1),
            Degrees::new(35.0),
            Degrees::new(80.0),
            Degrees::new(10.0),
            jd,
        ),
    );
    let demo_satellite_target = Target::new_static(demo_satellite.orbit.kepler_position(jd), jd);
    println!(
        "  {} target at JD {:.1}: r = {:.6}\n",
        demo_satellite.name,
        demo_satellite_target.time,
        demo_satellite_target.position.distance()
    );
}

fn section_target_with_proper_motion(jd: JulianDate) {
    println!("3) Target with proper motion (stellar-style target)");

    type MasPerYear = qtty::Per<qtty::MilliArcsecond, qtty::Year>;
    type MasPerYearQ = qtty::Quantity<MasPerYear>;

    let pm = ProperMotion::from_mu_alpha_star::<MasPerYear>(
        MasPerYearQ::new(27.54), // µα⋆
        MasPerYearQ::new(10.86), // µδ
    );

    let mut moving_target = Target::new(
        *catalog::BETELGEUSE.coordinate.get_position(),
        jd,
        pm.clone(),
    );

    println!(
        "  Betelgeuse-like target at J2000: RA {:.6}, Dec {:.6}",
        moving_target.position.ra(),
        moving_target.position.dec()
    );

    let jd_future = jd + 25.0 * JulianDate::JULIAN_YEAR;
    let moved = set_proper_motion_since_j2000(moving_target.position, pm, jd_future)
        .expect("proper-motion propagation should succeed");
    moving_target.update(moved, jd_future);

    println!(
        "  After 25 years: RA {:.6}, Dec {:.6}\n",
        moving_target.position.ra(),
        moving_target.position.dec()
    );
}

fn section_target_transform(jd: JulianDate) {
    println!("4) Target conversion across frame + center");

    type HelioEclAu = cartesian::position::EclipticMeanJ2000<AstronomicalUnit, Heliocentric>;
    type GeoEqAu = cartesian::position::EquatorialMeanJ2000<AstronomicalUnit, Geocentric>;

    let mars_helio: Target<HelioEclAu> = Target::new_static(Mars::vsop87a(jd), jd);
    let mars_geoeq: Target<GeoEqAu> = Target::from(&mars_helio);

    println!(
        "  Mars heliocentric ecliptic target: r = {:.6}",
        mars_helio.position.distance()
    );
    println!(
        "  Mars geocentric equatorial target: r = {:.6}",
        mars_geoeq.position.distance()
    );
}
