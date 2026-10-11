// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ESA LISA mission — OEM ephemerides and inter-satellite range
//!
//! This example demonstrates how to use reusable tabulated ephemerides and POD
//! observation helpers with LISA-like CCSDS OEM orbit files.
//!
//! ## What this example shows
//!
//! 1. Load one CCSDS OEM file per LISA spacecraft with
//!    [`siderust::ephemeris::tabulated::read_oem_tdb_ephemeris`].
//! 2. Build a typed [`siderust::ephemeris::tabulated::TabulatedEphemerisProvider`]
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

#![allow(clippy::print_stdout)]

use std::fs::File;
use std::path::Path;

use siderust::coordinates::centers::Heliocentric;
use siderust::coordinates::frames::EME2000;
use siderust::ephemeris::tabulated::{
    offset_tdb_epoch, read_oem_tdb_ephemeris, TabulatedEphemerisProvider,
};
use siderust::pod::observation::one_way_light_time_range;
use siderust::qtty::Second;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/test-data/lisa");
    let sc1 = read_oem_tdb_ephemeris(-1001, File::open(root.join("lisa_orbit_sample.oem1"))?)?;
    let sc2 = read_oem_tdb_ephemeris(-1002, File::open(root.join("lisa_orbit_sample.oem2"))?)?;
    let (start, _) = sc1.coverage()?;
    let provider = TabulatedEphemerisProvider::<Heliocentric, EME2000>::new(vec![sc1, sc2])?;

    let state = provider.state_at(-1001, start)?;
    println!("SC1 at coverage start: {}", state.position);
    println!("SC1 velocity: {}", state.velocity);

    let epoch = offset_tdb_epoch(start, Second::new(10.0))?;
    let range = one_way_light_time_range(&provider, -1001, -1002, epoch)?;
    println!("SC1 -> SC2 one-way range: {range:.3}");
    Ok(())
}
