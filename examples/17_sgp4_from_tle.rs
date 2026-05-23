//! Worked example: propagate a TLE through SGP4 to a few epochs.
//!
//! Loads the canonical Vallado SGP4 demo satellite (NORAD 5) via
//! [`siderust::formats::tle::parse_tle`], builds an [`siderust::astro::sgp4::Sgp4Propagator`], and prints the
//! TEME state at the TLE epoch and at three later epochs.
//!
//! Run with:
//! ```bash
//! cargo run --example 04_sgp4_from_tle
//! ```

#![allow(clippy::print_stdout, clippy::print_stderr)]

use qtty::Quantity;
use qtty::time::Day;
use siderust::formats::tle::{parse_tle, TleBuilder};
use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
use tempoch::{JulianDate, UTC};

const L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

fn main() {
    // 1. Parse a real TLE the conventional way (column-positional).
    let parsed = parse_tle(L1, L2).expect("vendored Vallado-VER TLE parses cleanly");

    // 2. Or re-synthesize the same record programmatically through the
    //    typed `TleBuilder` (proves the SGP4 pipeline accepts builder
    //    output, not just parser output).
    let synth = TleBuilder::new()
        .norad_id(parsed.norad_id)
        .international_designator(parsed.international_designator.clone())
        .classification(parsed.classification)
        .epoch(parsed.epoch)
        .inclination(parsed.inclination)
        .raan(parsed.raan)
        .eccentricity(parsed.eccentricity)
        .argument_of_perigee(parsed.argument_of_perigee)
        .mean_anomaly(parsed.mean_anomaly)
        .mean_motion(parsed.mean_motion)
        .mean_motion_dot(parsed.mean_motion_dot)
        .mean_motion_ddot(parsed.mean_motion_ddot)
        .bstar(parsed.bstar)
        .element_set_number(parsed.element_set_number)
        .revolution_number_at_epoch(parsed.revolution_number_at_epoch)
        .build()
        .expect("builder accepts round-tripped fields");

    // 3. Initialise the propagator (default WGS-72 gravity).
    let prop = Sgp4Propagator::from_tle_with_model(&synth, GravityModel::Wgs72)
        .expect("SGP4 init succeeds for NORAD 5");

    println!("Satellite     : NORAD {}", parsed.norad_id.0);
    println!("Gravity model : {:?}", prop.gravity_model());
    println!("TLE epoch (JD): {:.10}", prop.epoch_jd_utc().raw().value());
    println!();
    println!(
        "{:>12}  {:>14}  {:>14}  {:>14}  {:>10}  {:>10}  {:>10}",
        "t [min]", "x [km]", "y [km]", "z [km]", "vx [km/s]", "vy [km/s]", "vz [km/s]"
    );

    for &dt_min in &[0.0_f64, 360.0, 720.0, 1_080.0] {
        // Demonstrate both APIs: build a UTC Julian date, then call
        // `propagate_at`.
        let target = JulianDate::<UTC>::try_new(Quantity::<Day>::new(
            prop.epoch_jd_utc().raw().value() + dt_min / 1_440.0,
        ))
        .expect("finite Julian date");
        let s = prop.propagate_at(target).expect("propagation succeeds");
        let p = s.position().as_array();
        let v = s.velocity().as_array();
        println!(
            "{:>12.1}  {:>14.6}  {:>14.6}  {:>14.6}  {:>10.6}  {:>10.6}  {:>10.6}",
            dt_min,
            p[0].value(),
            p[1].value(),
            p[2].value(),
            v[0].value(),
            v[1].value(),
            v[2].value(),
        );
    }
}
