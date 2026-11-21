// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compare synthetic Gaia-style J2016.0 astrometry with a J2000.0 catalog.

#![allow(clippy::print_stdout)]

use siderust::catalogs::{parse_csv_chunk, CrossMatchOptions, LargeStarCatalog};
use siderust::qtty::{Arcsecond, Arcseconds, JulianYears};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut gaia = LargeStarCatalog::new();
    gaia.extend(parse_csv_chunk(
        "source_id,ra,dec,epoch,pmra,pmdec
1001,0.005,0,2016,1800,0
1002,32,10,2016,0,0
1003,80,-15,2016,0,0",
    )?);

    let mut historical = LargeStarCatalog::new();
    historical.extend(parse_csv_chunk(
        "source_id,ra,dec,epoch
1001,359.997,0,2000
2047,359.9971,0,2000
2048,32,10.00005,2000",
    )?);

    // By default, return all candidates and require complete motion for
    // off-epoch rows. Historical positions are already at J2000.0.
    let options = CrossMatchOptions::new(JulianYears::new(2000.0), Arcseconds::new(0.5).to());
    for candidate in gaia.cross_match(&historical, options)? {
        println!(
            "Gaia-style row {} (ID {}) → historical row {} (ID {}): {:.3} arcsec, {} candidates for Gaia row",
            candidate.left.index,
            candidate.left.source_id,
            candidate.right.index,
            candidate.right.source_id,
            candidate.separation.to::<Arcsecond>().value(),
            candidate.left_candidate_count,
        );
    }
    Ok(())
}
