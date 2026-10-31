// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Independent IMCCE check vectors rotated by ERFA, plus CSPICE DE440s states.
//! Regenerate with scripts/validate_ephemeris_frames.py; no kernel is bundled.
use serde_json::Value;
use siderust::bodies::solar_system::{Earth, Mars, Mercury, Moon, Neptune};
use siderust::coordinates::{
    cartesian::Position, centers::*, frames::*, transform::PositionAstroExt,
};
use siderust::ephemeris::{
    DynEphemeris, Ephemeris, MajorPlanet, PlanetEphemerisError, RuntimeEphemeris, Vsop87Ephemeris,
};
use siderust::qtty::*;
use siderust::time::{JulianDate, TDB, TT};

fn fixtures() -> Value {
    serde_json::from_str(include_str!("data/ephemeris_frames.json")).unwrap()
}
fn array(v: &Value) -> [f64; 3] {
    core::array::from_fn(|i| v[i].as_f64().unwrap())
}
fn close(actual: [f64; 3], expected: [f64; 3], eps: f64) {
    for i in 0..3 {
        assert!(
            (actual[i] - expected[i]).abs() < eps,
            "axis {i}: {actual:?} != {expected:?} (tolerance {eps})"
        );
    }
}
fn raw<C: siderust::coordinates::centers::ReferenceCenter<Params = ()>, U: LengthUnit>(
    p: Position<C, EclipticMeanJ2000, U>,
) -> [f64; 3] {
    [p.x().value(), p.y().value(), p.z().value()]
}
fn tt(jd_tdb: f64) -> JulianDate {
    tempoch::JulianDate::<TDB>::new(jd_tdb).to_scale::<TT>()
}

#[test]
fn published_vsop_check_states_in_canonical_axes() {
    let fixtures = fixtures();
    for case in fixtures["vsop"]
        .as_array()
        .unwrap()
        .iter()
        .chain(fixtures["vsop_endpoints"].as_array().unwrap())
    {
        // At the long-span endpoints, independent scalar argument evaluation
        // versus SIMD/FMA can differ at ~1e-10 AU; allow 1e-9 AU roundoff.
        let tolerance = if case["jd_tdb"].as_f64().unwrap() < 1_000_000.
            || case["jd_tdb"].as_f64().unwrap() > 3_000_000.
        {
            1e-9
        } else {
            1e-10
        };
        let jd = tt(case["jd_tdb"].as_f64().unwrap());
        macro_rules! state {
            ($body:ident) => {{
                let (p, v) = if case["version"] == "A" {
                    $body::vsop87a_pos_vel(jd)
                } else {
                    // Centers differ, so compare component arrays after this branch.
                    let (p, v) = $body::vsop87e_pos_vel(jd);
                    close(raw(p), array(&case["position"]), tolerance);
                    close(
                        [v.x().value(), v.y().value(), v.z().value()],
                        array(&case["velocity"]),
                        tolerance,
                    );
                    continue;
                };
                // IMCCE check file prints 10 decimals: 5e-11 AU (or AU/day)
                // quantization, plus summation and time conversion roundoff.
                close(raw(p), array(&case["position"]), tolerance);
                close(
                    [v.x().value(), v.y().value(), v.z().value()],
                    array(&case["velocity"]),
                    tolerance,
                );
            }};
        }
        match case["body"].as_str().unwrap() {
            "MERCURY" => state!(Mercury),
            "EARTH" => state!(Earth),
            "MARS" => state!(Mars),
            "NEPTUNE" => state!(Neptune),
            _ => unreachable!(),
        }
    }
}

#[test]
fn published_elp_table_h_in_canonical_axes() {
    for case in fixtures()["elp"].as_array().unwrap() {
        let jd = tt(case["jd_tdb"].as_f64().unwrap());
        // Table H is rounded to 1e-5 km. Allow 0.1 m for argument/summation
        // roundoff and the finite precision historical constants.
        close(
            raw(Moon::get_geo_position::<Kilometer>(jd)),
            array(&case["position"]),
            1e-4,
        );
    }
}

#[test]
fn positions_and_velocities_use_the_same_static_axes() {
    for days in [-365250., -36525., 0., 9789., 36525., 365250.] {
        let jd = siderust::J2000 + Days::new(days);
        let step = Days::new(0.01);
        macro_rules! derivative {
            ($body:ident, $position:ident, $velocity:ident) => {{
                let a = $body::$position(jd - step);
                let b = $body::$position(jd + step);
                let v = $body::$velocity(jd);
                let fd = [
                    (b.x() - a.x()).value() / 0.02,
                    (b.y() - a.y()).value() / 0.02,
                    (b.z() - a.z()).value() / 0.02,
                ];
                close(fd, [v.x().value(), v.y().value(), v.z().value()], 2e-8);
            }};
        }
        derivative!(Mercury, vsop87a, vsop87a_vel);
        derivative!(Earth, vsop87e, vsop87e_vel);
        derivative!(Neptune, vsop87a, vsop87a_vel);
    }
}

#[test]
fn normalized_position_roundtrip_preserves_vectors() {
    for days in [-365250., 0., 365250.] {
        let jd = siderust::J2000 + Days::new(days);
        let p = Earth::vsop87a(jd);
        let icrs = p.to_frame::<ICRS>(&jd);
        close(raw(icrs.to_frame::<EclipticMeanJ2000>(&jd)), raw(p), 1e-14);
    }
}

#[test]
#[ignore = "requires SIDERUST_BSP_PATH pointing to the official de440s.bsp"]
fn runtime_de440s_matches_cspice_erfa_reference() {
    let path = std::env::var("SIDERUST_BSP_PATH").expect("set SIDERUST_BSP_PATH");
    let eph = RuntimeEphemeris::from_bsp(path).unwrap();
    for case in fixtures()["jpl"].as_array().unwrap() {
        let jd = JulianDate::new(case["jd_tt"].as_f64().unwrap());
        let expected = array(&case["position"]);
        match (
            case["body"].as_i64().unwrap(),
            case["center"].as_i64().unwrap(),
        ) {
            (10, 0) => close(
                raw(eph.sun_barycentric(jd).to_unit::<Kilometer>()),
                expected,
                1e-6,
            ),
            (399, 0) => {
                close(
                    raw(eph.earth_barycentric(jd).to_unit::<Kilometer>()),
                    expected,
                    1e-6,
                );
                let v = eph
                    .earth_barycentric_velocity(jd)
                    .to_unit::<Per<Kilometer, Day>>();
                close(
                    [v.x().value(), v.y().value(), v.z().value()],
                    array(&case["velocity"]),
                    1e-7,
                );
            }
            (399, 10) => close(
                raw(eph.earth_heliocentric(jd).to_unit::<Kilometer>()),
                expected,
                1e-6,
            ),
            (301, 399) => close(raw(eph.moon_geocentric(jd)), expected, 1e-6),
            _ => unreachable!(),
        }
    }
}

// Deliberately unrelated planet state: proves center providers dispatch through
// the selected backend rather than consulting any built-in VSOP method.
struct TestEphemeris;
impl Ephemeris for TestEphemeris {
    fn sun_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        <Vsop87Ephemeris as Ephemeris>::sun_barycentric(jd)
    }
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        <Vsop87Ephemeris as Ephemeris>::earth_barycentric(jd)
    }
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        <Vsop87Ephemeris as Ephemeris>::earth_heliocentric(jd)
    }
    fn earth_barycentric_velocity(
        jd: JulianDate,
    ) -> siderust::coordinates::cartesian::Velocity<EclipticMeanJ2000, siderust::ephemeris::AuPerDay>
    {
        <Vsop87Ephemeris as Ephemeris>::earth_barycentric_velocity(jd)
    }
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        <Vsop87Ephemeris as Ephemeris>::moon_geocentric(jd)
    }
    fn try_planet_barycentric(
        _: MajorPlanet,
        _: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, PlanetEphemerisError>
    {
        Ok(Position::new(123., 456., 789.))
    }
}

#[test]
fn selected_backend_supplies_planetary_center_shifts() {
    use siderust::astro::{eop::NullEop, nutation::Iau2000B};
    use siderust::coordinates::transform::{AstroContext, CenterShiftProvider};
    let ctx = AstroContext::<TestEphemeris, NullEop>::default();
    let shift = <() as CenterShiftProvider<Marscentric, Barycentric, EclipticMeanJ2000>>::shift::<
        TestEphemeris,
        NullEop,
        Iau2000B,
    >(siderust::J2000, &ctx);
    close(shift.map(|x| x.value()), [123., 456., 789.], 1e-15);
    let inverse = <() as CenterShiftProvider<Barycentric, Marscentric, EclipticMeanJ2000>>::shift::<
        TestEphemeris,
        NullEop,
        Iau2000B,
    >(siderust::J2000, &ctx);
    close(inverse.map(|x| x.value()), [-123., -456., -789.], 1e-15);
    let geo = <() as CenterShiftProvider<Marscentric, Geocentric, EclipticMeanJ2000>>::shift::<
        TestEphemeris,
        NullEop,
        Iau2000B,
    >(siderust::J2000, &ctx);
    let earth = <TestEphemeris as Ephemeris>::earth_barycentric(siderust::J2000);
    close(
        geo.map(|x| x.value()),
        [
            123. - earth.x().value(),
            456. - earth.y().value(),
            789. - earth.z().value(),
        ],
        1e-13,
    );
}

#[test]
#[should_panic(expected = "selected ephemeris does not provide Pluto")]
fn missing_pluto_capability_cannot_silently_use_the_catalogue() {
    use siderust::astro::{eop::NullEop, nutation::Iau2000B};
    use siderust::coordinates::transform::{AstroContext, CenterShiftProvider};
    let ctx = AstroContext::<TestEphemeris, NullEop>::default();
    <() as CenterShiftProvider<Plutocentric, Barycentric, EclipticMeanJ2000>>::shift::<
        TestEphemeris,
        NullEop,
        Iau2000B,
    >(siderust::J2000, &ctx);
}

#[test]
fn default_pluto_center_keeps_its_documented_catalogue_model() {
    let jd = siderust::J2000;
    let p = <Vsop87Ephemeris as Ephemeris>::pluto_barycentric(jd);
    let h = siderust::bodies::solar_system::PLUTO
        .orbit
        .kepler_position(jd);
    let s = <Vsop87Ephemeris as Ephemeris>::sun_barycentric(jd);
    close(
        raw(p),
        [
            h.x().value() + s.x().value(),
            h.y().value() + s.y().value(),
            h.z().value() + s.z().value(),
        ],
        1e-14,
    );
}
