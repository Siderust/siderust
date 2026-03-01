// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime Ephemeris Example
//!
//! Demonstrates how to load a JPL DE4xx BSP file at runtime and query
//! planetary positions using [`RuntimeEphemeris`].
//!
//! Unlike the compile-time backends (DE440 via the `de440` feature),
//! DE441 and other large datasets are only supported at runtime because
//! embedding ~1.65 GB into the binary is impractical.
//!
//! ## Usage
//!
//! ```bash
//! # Load from a local BSP file
//! cargo run --example 12_runtime_ephemeris -- /path/to/de440.bsp
//!
//! # With the `runtime-data` feature: uses DataManager (auto-cache)
//! cargo run --features runtime-data --example 12_runtime_ephemeris -- /path/to/de440.bsp
//!
//! # With `runtime-data` and no path: shows how to download DE440 on first run
//! cargo run --features runtime-data --example 12_runtime_ephemeris
//! ```

use siderust::calculus::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::time::JulianDate;

fn print_positions(eph: &RuntimeEphemeris, label: &str) {
    let jd = JulianDate::J2000;

    let sun = eph.sun_barycentric(jd);
    let earth_bary = eph.earth_barycentric(jd);
    let earth_helio = eph.earth_heliocentric(jd);
    let earth_vel = eph.earth_barycentric_velocity(jd);
    let moon = eph.moon_geocentric(jd);

    println!("=== {} — positions at J2000 ===", label);
    println!("  Sun : {:.6}", sun);
    println!("  Earth : {:.6}", earth_bary);
    println!("  Earth : {:.6}", earth_helio);
    println!("  Earth vel: {:.8}", earth_vel);
    println!("  Moon : {:.1}", moon);
    println!();
}

/// Section 1: Load from a file path provided on the command line.
fn demo_load_from_path(path: &str) {
    println!("────────────────────────────────────────────────────");
    println!("1) Load from file: {}", path);
    println!("────────────────────────────────────────────────────");

    match RuntimeEphemeris::from_bsp(path) {
        Ok(eph) => {
            println!("  ✓ Loaded {:?}", eph);
            print_positions(&eph, "BSP file");
        }
        Err(e) => {
            eprintln!("  ✗ Failed to load '{}': {}", path, e);
            eprintln!("    Provide a valid DE440/DE441 BSP path as the first argument.");
        }
    }
}

/// Section 2: Load from bytes in memory (useful when embedding or downloading
/// the BSP data yourself before passing it to siderust).
fn demo_load_from_bytes() {
    println!("────────────────────────────────────────────────────");
    println!("2) Load from bytes in memory");
    println!("────────────────────────────────────────────────────");

    // In a real application you might read the file into a Vec<u8> and pass
    // it here — useful when you have already downloaded the BSP into a buffer.
    let fake_data = b"this is not a valid BSP file";
    match RuntimeEphemeris::from_bytes(fake_data) {
        Ok(_) => println!("  ✓ Parsed successfully"),
        Err(e) => println!("  ✗ Expected error for invalid data: {}", e),
    }
    println!();
}

/// Section 3 (requires `runtime-data` feature): Demonstrate [`DataManager`].
///
/// [`DataManager`] handles the persistent cache, integrity checking, and
/// (on explicit request) HTTP downloads from JPL servers.
#[cfg(feature = "runtime-data")]
fn demo_data_manager(explicit_download: bool) {
    use siderust::data::{DataManager, DatasetId};

    println!("────────────────────────────────────────────────────");
    println!("3) DataManager (feature = runtime-data)");
    println!("────────────────────────────────────────────────────");

    let dm = match DataManager::new() {
        Ok(dm) => {
            println!("  Cache dir: {}", dm.data_dir().display());
            dm
        }
        Err(e) => {
            eprintln!("  ✗ Cannot create DataManager: {}", e);
            return;
        }
    };

    // List all datasets and their cache status.
    println!("\n  Dataset availability:");
    for (id, cached) in dm.list() {
        let status = if cached {
            "✓ cached"
        } else {
            "✗ not cached"
        };
        println!("    {:10} — {}", format!("{:?}", id), status);
    }

    // DE440 (~120 MB): check if already downloaded before loading.
    let id = DatasetId::De440;
    println!("\n  Checking DE440...");

    if dm.is_available(id) {
        // File already cached — load it directly.
        match RuntimeEphemeris::from_data_manager(&dm, id) {
            Ok(eph) => print_positions(&eph, "DE440 (DataManager)"),
            Err(e) => eprintln!("  ✗ Load failed: {}", e),
        }
    } else if explicit_download {
        // User explicitly requested a download (--download flag).
        println!("  Downloading DE440 (~120 MB) from JPL...");
        println!("  NOTE: This is the ONLY way siderust downloads data —");
        println!("        it never happens automatically without your code calling it.");

        match dm.download(
            id,
            Some(Box::new(|downloaded, total| {
                if total > 0 {
                    let pct = downloaded * 100 / total;
                    eprint!(
                        "\r  Progress: {}% ({} / {} MB)",
                        pct,
                        downloaded >> 20,
                        total >> 20
                    );
                } else {
                    eprint!("\r  Downloaded: {} MB", downloaded >> 20);
                }
            })),
        ) {
            Ok(path) => {
                eprintln!(); // newline after progress
                println!("  ✓ Downloaded to: {}", path.display());
                match RuntimeEphemeris::from_bsp(&path) {
                    Ok(eph) => print_positions(&eph, "DE440 (freshly downloaded)"),
                    Err(e) => eprintln!("  ✗ Load failed: {}", e),
                }
            }
            Err(e) => eprintln!("\n  ✗ Download failed: {}", e),
        }
    } else {
        // Inform but do NOT download automatically — explicit consent required.
        println!("  DE440 is not cached yet.");
        println!();
        println!("  To download it, either:");
        println!("    a) Run this example with --download:         cargo run --features runtime-data --example 12_runtime_ephemeris -- --download");
        println!("    b) Call dm.ensure(DatasetId::De440) in code: it downloads if missing.");
        println!(
            "    c) Place de440.bsp in: {}",
            dm.data_dir().join("de440.bsp").display()
        );
        println!();
        println!("  DE441 (~1.65 GB) is only supported via runtime loading:");
        println!("    dm.ensure(DatasetId::De441)   // downloads de441_part-2.bsp on first call");
    }
    println!();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    println!("╔══════════════════════════════════════════════════╗");
    println!("║     Siderust Runtime Ephemeris Example           ║");
    println!("╚══════════════════════════════════════════════════╝");
    println!();

    // Demo 1: load from a user-provided file path (or show an error message).
    let bsp_path = args
        .iter()
        .skip(1)
        .find(|a| !a.starts_with("--"))
        .cloned()
        .unwrap_or_else(|| "de440.bsp".to_string());

    demo_load_from_path(&bsp_path);

    // Demo 2: loading from raw bytes.
    demo_load_from_bytes();

    // Demo 3: DataManager — only compiled when `runtime-data` feature is active.
    #[cfg(feature = "runtime-data")]
    {
        let explicit_download = args.iter().any(|a| a == "--download");
        demo_data_manager(explicit_download);
    }

    #[cfg(not(feature = "runtime-data"))]
    {
        println!("────────────────────────────────────────────────────");
        println!("3) DataManager (feature = runtime-data)  [SKIPPED]");
        println!("────────────────────────────────────────────────────");
        println!("  Enable with: cargo run --features runtime-data --example 12_runtime_ephemeris");
        println!();
    }
}
