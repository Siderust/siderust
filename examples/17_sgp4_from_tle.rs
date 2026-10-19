//! Propagate the Vallado demo satellite (NORAD 5) from its TLE using native Vallado SGP4.
//!
//! Run with: `cargo run --example 17_sgp4_from_tle`

#![allow(clippy::print_stdout)]

use siderust::astro::sgp4::Sgp4Propagator;
use siderust::formats::tle::parse_tle;
use siderust::qtty::Minutes;

const L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tle = parse_tle(L1, L2)?;
    let propagator = Sgp4Propagator::from_tle(&tle)?;

    println!("Satellite: NORAD {}", tle.norad_id.0);
    println!("Gravity model: {:?}", propagator.gravity_model());
    println!(
        "TLE epoch (UTC JD): {:.10}",
        propagator
            .epoch_jd_utc()
            .try_raw_with(propagator.time_context())?
    );

    for minutes in [0.0, 360.0, 720.0, 1_080.0] {
        let elapsed = Minutes::new(minutes);
        let state = propagator.propagate_minutes(elapsed)?;
        println!("\nElapsed: {elapsed:.1}");
        println!("TEME position: {:.6}", state.position());
        println!("TEME velocity: {:.6}", state.velocity());
    }
    Ok(())
}
