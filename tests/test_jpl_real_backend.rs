// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration smoke test for `RuntimeEphemeris` backed by a real BSP file.
//!
//! This test requires an actual JPL BSP file on disk.  It is skipped when the
//! `SIDERUST_BSP_PATH` environment variable is not set.
//!
//! Example:
//! ```bash
//! SIDERUST_BSP_PATH=/path/to/de440.bsp cargo test --test test_jpl_real_backend
//! ```

#[test]
fn runtime_ephemeris_real_bsp_smoke() {
    let Ok(path) = std::env::var("SIDERUST_BSP_PATH") else {
        eprintln!("Skipping RuntimeEphemeris smoke test: set SIDERUST_BSP_PATH to a BSP file.");
        return;
    };
    use siderust::calculus::ephemeris::{DynEphemeris, RuntimeEphemeris};
    use siderust::time::JulianDate;

    let eph = RuntimeEphemeris::from_bsp(&path)
        .unwrap_or_else(|e| panic!("Failed to load BSP file '{path}': {e}"));

    let earth = eph.earth_barycentric(JulianDate::J2000);
    assert!(
        earth.x().value().is_finite(),
        "RuntimeEphemeris: Earth barycentric should produce finite X at J2000"
    );
}
