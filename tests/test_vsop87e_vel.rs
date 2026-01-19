// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use siderust::astro::JulianDate;
use siderust::bodies::solar_system::{
    Earth, Jupiter, Mars, Mercury, Neptune, Saturn, Sun, Uranus, Venus,
};

const PRECISION: f64 = 1.0e-6;

macro_rules! check_body {
    ($body:ident) => {{
        let jd = JulianDate::J2000;
        let pos = $body::vsop87e(jd);
        let vel = $body::vsop87e_vel(jd);
        let (pos2, vel2) = $body::vsop87e_pos_vel(jd);
        let p1 = pos.get_position();
        let p2 = pos2.get_position();
        assert!((p1.x().value() - p2.x().value()).abs() < PRECISION);
        assert!((p1.y().value() - p2.y().value()).abs() < PRECISION);
        assert!((p1.z().value() - p2.z().value()).abs() < PRECISION);
        assert!((vel.x().value() - vel2.x().value()).abs() < PRECISION);
        assert!((vel.y().value() - vel2.y().value()).abs() < PRECISION);
        assert!((vel.z().value() - vel2.z().value()).abs() < PRECISION);
    }};
}

#[test]
fn velocities_match_combined() {
    check_body!(Mercury);
    check_body!(Venus);
    check_body!(Earth);
    check_body!(Mars);
    check_body!(Jupiter);
    check_body!(Saturn);
    check_body!(Uranus);
    check_body!(Neptune);
}

#[test]
fn sun_position_finite() {
    let jd = JulianDate::J2000;
    let pos = *Sun::vsop87e(jd).get_position();
    assert!(pos.x().value().is_finite());
    assert!(pos.y().value().is_finite());
    assert!(pos.z().value().is_finite());
}
