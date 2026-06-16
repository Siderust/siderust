// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ESA LISA mission — OEM ephemerides and inter-satellite range
//!
//! This example demonstrates how to use the reusable `siderust::pod` tabulated
//! ephemeris and observation helpers with LISA-like CCSDS OEM orbit files.
//!
//! ## What this example shows
//!
//! 1. Load one CCSDS OEM file per LISA spacecraft with
//!    [`siderust::pod::tabulated::read_oem_tdb_ephemeris`].
//! 2. Build a typed [`siderust::pod::tabulated::TabulatedEphemerisProvider`]
//!    backed by cubic Hermite interpolation.
//! 3. Query a typed spacecraft state and compute a one-way light-time-corrected
//!    inter-spacecraft range with
//!    [`siderust::pod::observation::one_way_light_time_range`].
//!
//! ## Running the example
//!
//! ```text
//! cargo run --example 18_lisa_pod --features pod
//! ```
//!
//! The example loads orbit files from `tests/test-data/lisa/`.
//!
//! ## References
//!
//! - Consultative Committee for Space Data Systems. (2019). Orbit Data
//!   Messages, CCSDS 502.0-B-3.
//! - Martens, W., Joffre, E. (2021). Trajectory Design for the ESA LISA
//!   Mission. *J. Astronaut. Sci.*, 68, 402-443.
//!   <https://doi.org/10.1007/s40295-021-00263-2>
//! - Danzmann, K. et al. (2017). *LISA: Laser Interferometer Space Antenna.*
//!   ESA/SRE(2017)1. <https://arxiv.org/abs/1702.00786>

#![allow(clippy::print_stdout, missing_docs)]

use std::fs::File;

use qtty::Second;
use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::EME2000;
use siderust::pod::observation::one_way_light_time_range;
use siderust::pod::tabulated::{
    read_oem_tdb_ephemeris, TabulatedEphemeris, TabulatedEphemerisProvider,
};
use tempoch::{Time, TDB};

type LisaEphemeris = TabulatedEphemeris<Heliocentric, EME2000>;
type LisaProvider = TabulatedEphemerisProvider<Heliocentric, EME2000>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let provider = lisa_provider()?;
    let (coverage_start, _) = provider
        .ephemeris(-1001)
        .expect("SC1 ephemeris")
        .coverage()?;

    let sc1 = provider.state_at(-1001, coverage_start)?;
    println!(
        "SC1 at coverage start: ({:.0}, {:.0}, {:.0}) km",
        sc1.position.x().value(),
        sc1.position.y().value(),
        sc1.position.z().value()
    );
    println!(
        "SC1 velocity: ({:.1}, {:.1}, {:.1}) km/s",
        sc1.velocity.x().value(),
        sc1.velocity.y().value(),
        sc1.velocity.z().value()
    );

    let range_epoch = Time::<TDB>::from_raw_j2000_seconds(Second::new(
        coverage_start.to::<tempoch::J2000s>().raw().value() + 10.0,
    ))?;
    let measured = one_way_light_time_range(&provider, -1001, -1002, range_epoch)?;
    let modeled = one_way_light_time_range(&provider, -1001, -1002, range_epoch)?;
    let residual_m = measured.value() - modeled.value();

    println!("SC1 -> SC2 one-way range: {:.3} m", modeled.value());
    println!("synthetic O-C residual: {residual_m:.6} m");

    assert!(residual_m.abs() < 1e-9);
    Ok(())
}

fn lisa_provider() -> Result<LisaProvider, Box<dyn std::error::Error>> {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test-data/lisa");
    Ok(LisaProvider::new(vec![
        load_ephemeris(-1001, &format!("{root}/lisa_orbit_sample.oem1"))?,
        load_ephemeris(-1002, &format!("{root}/lisa_orbit_sample.oem2"))?,
        load_ephemeris(-1003, &format!("{root}/lisa_orbit_sample.oem3"))?,
    ])?)
}

fn load_ephemeris(
    body_naif_id: i32,
    path: &str,
) -> Result<LisaEphemeris, Box<dyn std::error::Error>> {
    Ok(read_oem_tdb_ephemeris(body_naif_id, File::open(path)?)?)
}
