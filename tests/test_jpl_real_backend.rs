// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integration guard for real (non-stubbed) JPL backends.

#[cfg(all(feature = "de441", not(siderust_mock_de441)))]
#[test]
fn de441_real_backend_smoke() {
    use siderust::calculus::ephemeris::{De441Ephemeris, Ephemeris};
    use siderust::time::JulianDate;

    let backend_type = core::any::type_name::<De441Ephemeris>();
    assert!(
        backend_type.contains("DeEphemeris"),
        "Expected real DE441 backend, got {backend_type}"
    );

    let earth = De441Ephemeris::earth_barycentric(JulianDate::J2000);
    assert!(
        earth.x().value().is_finite(),
        "DE441 real backend should produce finite coordinates at J2000"
    );
}
