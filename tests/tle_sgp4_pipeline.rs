// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! End-to-end smoke test for the typed TLE → SGP4 → TEME pipeline.
//!
//! This is the canonical "the typed boundary holds" check for
//! `siderust-pod`. It exercises:
//!
//! 1. [`parse_tle`] consuming a published Vallado SGP4 Verification TLE.
//! 2. [`Sgp4Propagator::from_tle_with_model`] adapting the TLE to the
//!    Siderust-native propagation core with the WGS-84 gravity model.
//! 3. Propagation at the TLE epoch and at later epochs, returning
//!    typed [`TemeState`] values whose components live in the typed
//!    [`affn`] / [`qtty`] / [`tempoch`] coordinate stack.
//!
//! No external dataset is required. If this test fails the typed
//! boundary between `tle`, `sgp4`, and the `siderust` / `tempoch` /
//! `qtty` stack is broken.
//!
//! [`parse_tle`]: siderust::formats::tle::parse_tle
//! [`Sgp4Propagator::from_tle_with_model`]: siderust::astro::sgp4::Sgp4Propagator::from_tle_with_model
//! [`TemeState`]: siderust::astro::sgp4::TemeState

use qtty::time::Day;
use qtty::Quantity;
use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
use siderust::formats::tle::parse_tle;
use tempoch::JulianDate;

// Vallado SGP4 Verification TLE, NORAD catalog 5.
const L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

#[test]
fn tle_sgp4_pipeline_produces_typed_teme_state() {
    let tle = parse_tle(L1, L2).expect("vendored Vallado-VER TLE parses cleanly");
    let prop = Sgp4Propagator::from_tle_with_model(&tle, GravityModel::Wgs84)
        .expect("Vallado-VER TLE accepted by the native propagator");

    let epoch_jd: JulianDate<_> = tle.epoch.to::<tempoch::JD>();
    let state_at_epoch = prop
        .propagate_at(epoch_jd)
        .expect("propagation at the TLE epoch must succeed");

    // The TEME position lives in the typed affn/qtty stack: components
    // are typed quantities in kilometres, and the radial magnitude
    // for a LEO test satellite must be a physically plausible Earth-
    // centred distance (≥ 6 000 km).
    let x = state_at_epoch.position().x().value();
    let y = state_at_epoch.position().y().value();
    let z = state_at_epoch.position().z().value();
    let radius_km = (x * x + y * y + z * z).sqrt();
    assert!(
        (6_000.0..200_000.0).contains(&radius_km),
        "TEME radius {radius_km} km is not a plausible Earth-centred distance"
    );

    // Velocity is finite and non-zero.
    let vx = state_at_epoch.velocity().x().value();
    let vy = state_at_epoch.velocity().y().value();
    let vz = state_at_epoch.velocity().z().value();
    let speed_kms = (vx * vx + vy * vy + vz * vz).sqrt();
    assert!(
        (0.5..15.0).contains(&speed_kms),
        "TEME speed {speed_kms} km/s outside any physical envelope for an Earth satellite"
    );
}

#[test]
fn tle_sgp4_pipeline_propagates_forward_in_time() {
    let tle = parse_tle(L1, L2).expect("vendored Vallado-VER TLE parses cleanly");
    let prop =
        Sgp4Propagator::from_tle_with_model(&tle, GravityModel::Wgs84).expect("propagator builds");

    let epoch_jd: JulianDate<_> = tle.epoch.to::<tempoch::JD>();
    let state_0 = prop.propagate_at(epoch_jd).expect("epoch state");
    // Advance by one minute (TLE epoch is a JulianDate; adding a typed
    // Day quantity exercises the typed time addition path).
    let one_minute = Quantity::<Day>::new(1.0 / 1440.0);
    let later_jd = epoch_jd + one_minute;
    let state_1 = prop.propagate_at(later_jd).expect("later state");

    let dx = state_1.position().x().value() - state_0.position().x().value();
    let dy = state_1.position().y().value() - state_0.position().y().value();
    let dz = state_1.position().z().value() - state_0.position().z().value();
    let moved_km = (dx * dx + dy * dy + dz * dz).sqrt();
    // At LEO orbital speeds the spacecraft sweeps at least a few
    // hundred kilometres in one minute.
    assert!(
        moved_km > 50.0,
        "spacecraft barely moved ({moved_km} km) in one minute; pipeline likely stuck at epoch"
    );
}
