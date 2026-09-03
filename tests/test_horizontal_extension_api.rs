// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Public-API regression tests for the ergonomic horizontal-coordinate extensions.
//!
//! These tests intentionally exercise the convenience methods added by PR #86,
//! rather than calling the lower-level horizontal transformation traits directly.

use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Topocentric;
use siderust::coordinates::frames::EquatorialTrueOfDate;
use siderust::coordinates::spherical;
use siderust::coordinates::transform::{
    DirectionAstroExt, SphericalDirectionAstroExt, TopocentricEquatorialExt,
    TopocentricHorizontalExt,
};
use siderust::qtty::*;
use siderust::time::JulianDate;

const ANGULAR_EPSILON_DEG: f64 = 1e-9;
const POSITION_EPSILON: f64 = 1e-10;

#[test]
fn cartesian_direction_horizontal_extensions_round_trip() {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let jd_ut1 = JulianDate::new(2_459_015.5);
    let jd_tt = jd_ut1.add_seconds(69.184);
    let equatorial = spherical::direction::EquatorialTrueOfDate::new(
        Degrees::new(101.287),
        Degrees::new(-16.716),
    )
    .to_cartesian();

    // TT-only convenience path: UT1 is derived from the default EOP context.
    let horizontal = DirectionAstroExt::to_horizontal(&equatorial, &jd_tt, &site);
    let back = DirectionAstroExt::to_equatorial(&horizontal, &jd_tt, &site);
    let separation = equatorial
        .to_spherical()
        .angular_separation(&back.to_spherical())
        .value();
    assert!(
        separation < ANGULAR_EPSILON_DEG,
        "cartesian TT-only round-trip separation: {separation} deg"
    );

    // Explicit UT1+TT path.
    let horizontal = DirectionAstroExt::to_horizontal_precise(&equatorial, &jd_tt, &jd_ut1, &site);
    let back = DirectionAstroExt::to_equatorial_precise(&horizontal, &jd_tt, &jd_ut1, &site);
    let separation = equatorial
        .to_spherical()
        .angular_separation(&back.to_spherical())
        .value();
    assert!(
        separation < ANGULAR_EPSILON_DEG,
        "cartesian precise round-trip separation: {separation} deg"
    );
}

#[test]
fn spherical_direction_horizontal_extensions_round_trip() {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let jd_ut1 = JulianDate::new(2_459_015.5);
    let jd_tt = jd_ut1.add_seconds(69.184);
    let equatorial = spherical::direction::EquatorialTrueOfDate::new(
        Degrees::new(101.287),
        Degrees::new(-16.716),
    );

    // TT-only spherical wrapper.
    let horizontal = SphericalDirectionAstroExt::to_horizontal(&equatorial, &jd_tt, &site);
    let back = SphericalDirectionAstroExt::to_equatorial(&horizontal, &jd_tt, &site);
    let separation = equatorial.angular_separation(&back).value();
    assert!(
        separation < ANGULAR_EPSILON_DEG,
        "spherical TT-only round-trip separation: {separation} deg"
    );

    // Explicit UT1+TT spherical wrapper.
    let horizontal =
        SphericalDirectionAstroExt::to_horizontal_precise(&equatorial, &jd_tt, &jd_ut1, &site);
    let back =
        SphericalDirectionAstroExt::to_equatorial_precise(&horizontal, &jd_tt, &jd_ut1, &site);
    let separation = equatorial.angular_separation(&back).value();
    assert!(
        separation < ANGULAR_EPSILON_DEG,
        "spherical precise round-trip separation: {separation} deg"
    );
}

#[test]
fn topocentric_position_aliases_match_explicit_methods() {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let jd_ut1 = JulianDate::new(2_459_015.5);
    let jd_tt = jd_ut1.add_seconds(69.184);
    let distance = 1_000.0 * KM;

    let equatorial_direction = spherical::direction::EquatorialTrueOfDate::new(
        Degrees::new(101.287),
        Degrees::new(-16.716),
    );
    let equatorial_spherical =
        equatorial_direction.position_with_params::<Topocentric, Kilometer>(site, distance);
    let equatorial = Position::<Topocentric, EquatorialTrueOfDate, Kilometer>::from_spherical(
        &equatorial_spherical,
    );

    let horizontal_explicit =
        TopocentricEquatorialExt::to_horizontal_position(&equatorial, &jd_ut1, &jd_tt);
    let horizontal_alias = TopocentricEquatorialExt::to_horizontal(&equatorial, &jd_ut1, &jd_tt);

    assert_eq!(horizontal_alias.center_params(), &site);
    assert!(
        (horizontal_alias.x().value() - horizontal_explicit.x().value()).abs() < POSITION_EPSILON
    );
    assert!(
        (horizontal_alias.y().value() - horizontal_explicit.y().value()).abs() < POSITION_EPSILON
    );
    assert!(
        (horizontal_alias.z().value() - horizontal_explicit.z().value()).abs() < POSITION_EPSILON
    );

    let equatorial_explicit =
        TopocentricHorizontalExt::to_equatorial_position(&horizontal_alias, &jd_ut1, &jd_tt);
    let equatorial_alias =
        TopocentricHorizontalExt::to_equatorial(&horizontal_alias, &jd_ut1, &jd_tt);

    assert_eq!(equatorial_alias.center_params(), &site);
    assert!(
        (equatorial_alias.x().value() - equatorial_explicit.x().value()).abs() < POSITION_EPSILON
    );
    assert!(
        (equatorial_alias.y().value() - equatorial_explicit.y().value()).abs() < POSITION_EPSILON
    );
    assert!(
        (equatorial_alias.z().value() - equatorial_explicit.z().value()).abs() < POSITION_EPSILON
    );
}
