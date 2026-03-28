// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Nutation Model Selection Example
//!
//! Shows the two intended usage modes:
//! - default transforms with no explicit context
//! - custom transforms using `AstroContext::with_model::<...>()`
//!
//! Run with: `cargo run --example 14_nutation_models`

use siderust::astro::nutation::{Iau2000B, Iau2006};
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::frames::{EquatorialTrueOfDate, ICRS};
use siderust::coordinates::transform::{AstroContext, DirectionAstroExt};
use siderust::time::JulianDate;

fn chord_delta(a: &Direction<EquatorialTrueOfDate>, b: &Direction<EquatorialTrueOfDate>) -> f64 {
    ((a.x() - b.x()).powi(2) + (a.y() - b.y()).powi(2) + (a.z() - b.z()).powi(2)).sqrt()
}

fn main() {
    let jd = JulianDate::new(2_458_850.0);
    let icrs = Direction::<ICRS>::new(0.6, -0.3, 0.74);

    println!("=== Nutation Model Selection ===\n");
    println!("Epoch (TT): {jd}");
    println!(
        "Input ICRS direction: ({:.6}, {:.6}, {:.6})\n",
        icrs.x(),
        icrs.y(),
        icrs.z()
    );

    // Default path: no explicit context. This uses DefaultNutationModel = Iau2006A.
    let default_true: Direction<EquatorialTrueOfDate> = icrs.to_frame(&jd);

    // Custom path: bind a specific nutation model to the runtime context.
    let ctx = AstroContext::new();
    let fast_ctx = ctx.with_model::<Iau2000B>();
    let precession_only_ctx = ctx.with_model::<Iau2006>();

    let fast_true: Direction<EquatorialTrueOfDate> = icrs.to_frame_with(&jd, &fast_ctx);
    let precession_only_true: Direction<EquatorialTrueOfDate> =
        icrs.to_frame_with(&jd, &precession_only_ctx);

    println!("1. Default transform (Iau2006A)");
    println!(
        "   TOD direction = ({:.12}, {:.12}, {:.12})",
        default_true.x(),
        default_true.y(),
        default_true.z()
    );

    println!("\n2. Custom transform with Iau2000B");
    println!(
        "   TOD direction = ({:.12}, {:.12}, {:.12})",
        fast_true.x(),
        fast_true.y(),
        fast_true.z()
    );
    println!(
        "   Chord delta vs default = {:.3e}",
        chord_delta(&default_true, &fast_true)
    );

    println!("\n3. Custom transform with Iau2006 (precession only)");
    println!(
        "   TOD direction = ({:.12}, {:.12}, {:.12})",
        precession_only_true.x(),
        precession_only_true.y(),
        precession_only_true.z()
    );
    println!(
        "   Chord delta vs default = {:.3e}",
        chord_delta(&default_true, &precession_only_true)
    );

    println!("\nPattern to copy:");
    println!("  let ctx = AstroContext::new();");
    println!("  let custom = ctx.with_model::<Iau2000B>();");
    println!("  let tod: Direction<EquatorialTrueOfDate> = icrs.to_frame_with(&jd, &custom);");
}
