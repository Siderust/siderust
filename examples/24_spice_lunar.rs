// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon
//! Geometric DE440 Sun/Earth observations; see doc/spice_lunar.md for kernels.
#![allow(clippy::print_stdout)]
use siderust::qtty::{DEG, M};
use siderust::spice::{LunarObserver, LunarTarget, SpiceContext};
use siderust::tempoch::{JulianDate, TDB, TT};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let directory = std::env::args()
        .nth(1)
        .ok_or("Usage: 24_spice_lunar KERNEL_DIRECTORY")?;
    let mut kernels = SpiceContext::new();
    for file in [
        "de440s.bsp",
        "moon_pa_de440_200625.bpc",
        "moon_de440_220930.tf",
        "pck00011.tpc",
        "naif0012.tls",
    ] {
        kernels.load(std::path::Path::new(&directory).join(file))?;
    }
    // Actual radial distance: east-positive lon, planetocentric lat, radial up.
    // Site coordinates belong to MOON_ME_DE440_ME421 via this FK's MOON_ME.
    let observer = LunarObserver::radial(
        &kernels,
        "MOON_ME",
        32.3 * DEG,
        -89.0 * DEG,
        1_738_000.0 * M,
    )?;
    // Explicit TT -> TDB conversion. Direct CSPICE parity uses identical ET.
    let epoch = JulianDate::<TT>::new(2_461_358.5).to::<TDB>().to_j2000s();
    for target in [LunarTarget::Sun, LunarTarget::Earth] {
        let result = observer.observe(target, epoch)?;
        println!(
            "{target:?}: alt={}, az={}, range={}, above={}",
            result.horizontal.alt(),
            result.horizontal.az(),
            result.horizontal.range(),
            result.above_horizon()
        );
    }
    Ok(())
}
