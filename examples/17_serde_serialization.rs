// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Serde serialization examples.
//!
//! Run with: `cargo run --example 17_serde_serialization --features serde`

use qtty::*;
use serde::{Deserialize, Serialize};
use siderust::astro::orbit::Orbit;
use siderust::bodies::comet::HALLEY;
use siderust::bodies::solar_system::{Earth, Mars, Moon};
use siderust::coordinates::{cartesian, centers, frames, spherical};
use siderust::targets::{Target, Trackable};
use siderust::time::{JulianDate, ModifiedJulianDate};
use std::fs;

#[derive(Debug, Serialize, Deserialize)]
struct TimeBundle {
    j2000: JulianDate,
    mjd: ModifiedJulianDate,
    timeline: Vec<JulianDate>,
}

#[derive(Debug, Serialize, Deserialize)]
struct CoordinateBundle {
    geo_icrs_cart:
        cartesian::Position<centers::Geocentric, frames::ICRS, Kilometer>,
    helio_ecl_sph:
        spherical::Position<centers::Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>,
    observer_site: centers::Geodetic<frames::ECEF>,
}

#[derive(Debug, Serialize, Deserialize)]
struct BodySnapshot {
    name: String,
    epoch: JulianDate,
    orbit: Orbit,
    heliocentric_ecliptic:
        cartesian::Position<centers::Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>,
}

#[derive(Debug, Serialize, Deserialize)]
struct BodyTargetsBundle {
    // Mars `Trackable` output is already a CoordinateWithPM<...>, i.e. Target<...>.
    mars_bary_target:
        Target<cartesian::Position<centers::Barycentric, frames::EclipticMeanJ2000, AstronomicalUnit>>,
    // Moon does not return CoordinateWithPM from `track`, so we wrap a snapshot.
    moon_geo_target:
        Target<cartesian::Position<centers::Geocentric, frames::EclipticMeanJ2000, Kilometer>>,
}

fn pretty_json<T: Serialize>(value: &T) -> String {
    serde_json::to_string_pretty(value).expect("serialize to pretty JSON")
}

fn roundtrip<T>(value: &T) -> T
where
    T: Serialize + for<'de> Deserialize<'de>,
{
    let json = serde_json::to_string(value).expect("serialize");
    serde_json::from_str(&json).expect("deserialize")
}

fn main() {
    println!("=== Siderust Serde Serialization Examples ===\n");

    let jd = JulianDate::J2000;

    // =========================================================================
    // 1) Times
    // =========================================================================
    println!("1) TIME OBJECTS");
    println!("---------------");

    let time_bundle = TimeBundle {
        j2000: jd,
        mjd: ModifiedJulianDate::from(jd),
        timeline: vec![jd, jd + Days::new(1.0), jd + Days::new(7.0)],
    };
    println!("{}", pretty_json(&time_bundle));

    let recovered_times: TimeBundle = roundtrip(&time_bundle);
    println!(
        "Roundtrip check: j2000={:.1}, timeline_len={}\n",
        recovered_times.j2000.value(),
        recovered_times.timeline.len()
    );

    // =========================================================================
    // 2) Coordinates
    // =========================================================================
    println!("2) COORDINATE OBJECTS");
    println!("---------------------");

    let coords = CoordinateBundle {
        geo_icrs_cart:
            cartesian::Position::<centers::Geocentric, frames::ICRS, Kilometer>::new(
                6371.0, 0.0, 0.0,
            ),
        helio_ecl_sph:
            spherical::Position::<centers::Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>::new_raw(
                Degrees::new(5.0),   // lat
                Degrees::new(120.0), // lon
                AstronomicalUnits::new(1.2),
            ),
        observer_site:
            centers::Geodetic::<frames::ECEF>::new(
                Degrees::new(-17.8947), // lon
                Degrees::new(28.7636),  // lat
                Meters::new(2396.0),    // height
            ),
    };
    println!("{}", pretty_json(&coords));

    let recovered_coords: CoordinateBundle = roundtrip(&coords);
    println!(
        "Roundtrip check: x={:.1} km, lon={:.4} deg\n",
        recovered_coords.geo_icrs_cart.x().value(),
        recovered_coords.observer_site.lon.value()
    );

    // =========================================================================
    // 3) Body-related objects: orbit + ephemeris snapshots
    // =========================================================================
    println!("3) BODY-RELATED OBJECTS");
    println!("-----------------------");

    // NOTE:
    // `Planet`/`Star` structs are not serde-derived in the current API.
    // We serialize body-related data that *is* serde-ready: orbit elements
    // and concrete coordinate snapshots.
    let earth_snapshot = BodySnapshot {
        name: "Earth".to_string(),
        epoch: jd,
        orbit: siderust::bodies::EARTH.orbit,
        heliocentric_ecliptic: Earth::vsop87a(jd),
    };
    let halley_snapshot = BodySnapshot {
        name: HALLEY.name.to_string(),
        epoch: jd,
        orbit: HALLEY.orbit,
        heliocentric_ecliptic: HALLEY.orbit.kepler_position(jd),
    };

    println!("Earth snapshot JSON:");
    println!("{}", pretty_json(&earth_snapshot));
    println!("Halley snapshot JSON:");
    println!("{}", pretty_json(&halley_snapshot));

    let recovered_halley: BodySnapshot = roundtrip(&halley_snapshot);
    println!(
        "Roundtrip check: {} @ JD {:.1}, r={:.6} AU\n",
        recovered_halley.name,
        recovered_halley.epoch.value(),
        recovered_halley.heliocentric_ecliptic.distance().value()
    );

    // =========================================================================
    // 4) Target objects (CoordinateWithPM alias)
    // =========================================================================
    println!("4) TARGET OBJECTS");
    println!("-----------------");

    let mars_target = Mars.track(jd);
    let moon_target = Target::new_static(Moon.track(jd), jd);

    let targets = BodyTargetsBundle {
        mars_bary_target: mars_target,
        moon_geo_target: moon_target,
    };
    println!("{}", pretty_json(&targets));

    let recovered_targets: BodyTargetsBundle = roundtrip(&targets);
    println!(
        "Roundtrip check: Mars target JD {:.1}, Moon target JD {:.1}\n",
        recovered_targets.mars_bary_target.time.value(),
        recovered_targets.moon_geo_target.time.value()
    );

    // =========================================================================
    // 5) File I/O
    // =========================================================================
    println!("5) FILE I/O");
    println!("----------");

    let out_path = "/tmp/siderust_serde_example_targets.json";
    fs::write(out_path, pretty_json(&targets)).expect("write JSON file");
    let loaded = fs::read_to_string(out_path).expect("read JSON file");
    let _: BodyTargetsBundle = serde_json::from_str(&loaded).expect("deserialize loaded JSON");
    println!("Saved and loaded: {}", out_path);
}
