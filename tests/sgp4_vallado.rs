// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Independent published Vallado vectors, exercised through the typed public API.
//! Fixture provenance and license: tests/fixtures/sgp4/README.md.

use siderust::astro::sgp4::{GravityModel, Sgp4Error, Sgp4Propagator, TemeState};
use siderust::formats::tle::{parse_tle, TLE};
use siderust::qtty::{Degrees, Minutes};
use tempoch::{Time, TimeContext, JD, UTC};

const L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

fn assert_state(state: &TemeState, position: [f64; 3], velocity: [f64; 3], p: f64, v: f64) {
    for (i, (actual, expected)) in state.position().as_array().iter().zip(position).enumerate() {
        assert!(
            (actual.value() - expected).abs() < p,
            "position[{i}]: {} vs {expected}",
            actual.value()
        );
    }
    for (i, (actual, expected)) in state.velocity().as_array().iter().zip(velocity).enumerate() {
        assert!(
            (actual.value() - expected).abs() < v,
            "velocity[{i}]: {} vs {expected}",
            actual.value()
        );
    }
}

#[test]
fn full_vallado_verification_corpus() {
    let mut lines = include_str!("fixtures/sgp4/vallado.txt")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .peekable();
    let mut cases = 0;
    let mut states = 0;
    while let Some(line1) = lines.next() {
        let line2 = lines.next().unwrap();
        let tle = parse_tle(line1, line2).unwrap();
        let prop = Sgp4Propagator::from_tle(&tle).unwrap();
        cases += 1;
        while lines.peek().is_some_and(|l| !l.starts_with("1 ")) {
            let values: Vec<_> = lines.next().unwrap().split_whitespace().collect();
            let minutes = values[0].parse::<f64>().unwrap();
            let result = prop.propagate_minutes(Minutes::new(minutes));
            states += 1;
            if values[1] == "ERROR" {
                let error = result.unwrap_err();
                assert!(matches!(error, Sgp4Error::Propagation { .. }));
                assert!(error.to_string().contains(values[2]));
            } else {
                let numbers: Vec<f64> = values[1..].iter().map(|s| s.parse().unwrap()).collect();
                let state = result.unwrap();
                assert_state(
                    &state,
                    numbers[..3].try_into().unwrap(),
                    numbers[3..].try_into().unwrap(),
                    1e-6,
                    1e-9,
                );
                // The epoch produced by nominal-minute propagation retains the
                // instant in compensated Tempoch storage, avoiding scalar JD loss.
                let at = prop.propagate_at(state.epoch()).unwrap();
                assert_state(
                    &at,
                    state.position().as_array().map(|x| x.value()),
                    state.velocity().as_array().map(|x| x.value()),
                    1e-6,
                    1e-9,
                );
            }
        }
    }
    assert_eq!(cases, 33);
    assert_eq!(states, 668);
}

#[test]
fn rejects_invalid_mutated_elements() {
    let valid = parse_tle(L1, L2).unwrap();
    let mut invalid: Vec<TLE> = Vec::new();
    for e in [-0.1, 1.0, f64::NAN, f64::INFINITY] {
        let mut tle = valid.clone();
        tle.eccentricity = e;
        invalid.push(tle);
    }
    for n in [
        -1.0,
        0.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        let mut tle = valid.clone();
        tle.mean_motion = qtty::Quantity::new(n);
        invalid.push(tle);
    }
    for i in [-1.0, 181.0, f64::NAN] {
        let mut tle = valid.clone();
        tle.inclination = Degrees::new(i);
        invalid.push(tle);
    }
    for field in 0..4 {
        let mut tle = valid.clone();
        match field {
            0 => tle.raan = Degrees::new(f64::NAN),
            1 => tle.argument_of_perigee = Degrees::new(f64::INFINITY),
            2 => tle.mean_anomaly = Degrees::new(f64::NAN),
            _ => tle.bstar = f64::INFINITY,
        }
        invalid.push(tle);
    }
    for tle in invalid {
        assert!(matches!(
            Sgp4Propagator::from_tle(&tle),
            Err(Sgp4Error::InvalidElements { .. })
        ));
    }
}

#[test]
fn invalid_time_offsets_and_epochs_return_errors() {
    let mut tle = parse_tle(L1, L2).unwrap();
    let prop = Sgp4Propagator::from_tle(&tle).unwrap();
    for offset in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MAX,
        -f64::MAX,
    ] {
        assert!(matches!(
            prop.propagate_minutes(Minutes::new(offset)),
            Err(Sgp4Error::TimeConversion(_))
        ));
    }
    // Tempoch allows an infinite boundary instant, but civil UTC cannot encode it.
    let infinite = Time::<UTC>::from_raw_j2000_seconds(qtty::Second::new(f64::INFINITY)).unwrap();
    tle.epoch = infinite;
    assert!(matches!(
        Sgp4Propagator::from_tle(&tle),
        Err(Sgp4Error::TimeConversion(_))
    ));
    assert!(matches!(
        prop.propagate_at(infinite.reinterpret()),
        Err(Sgp4Error::TimeConversion(_))
    ));
}

#[test]
fn utc_nominal_minutes_cross_leap_seconds() {
    let ctx = TimeContext::new();
    let mut tle = parse_tle(L1, L2).unwrap();
    tle.epoch =
        Time::<UTC>::try_from_chrono_with("2016-12-31T23:59:00Z".parse().unwrap(), &ctx).unwrap();
    let prop = Sgp4Propagator::from_tle(&tle).unwrap();
    for (label, minutes) in [
        ("2017-01-01T00:00:00Z", 1.0),
        ("2016-12-30T23:59:00Z", -1440.0),
    ] {
        let target = Time::<UTC>::try_from_chrono_with(label.parse().unwrap(), &ctx)
            .unwrap()
            .to_with::<JD>(&ctx)
            .unwrap();
        let at = prop.propagate_at(target).unwrap();
        let offset = prop.propagate_minutes(Minutes::new(minutes)).unwrap();
        assert_eq!(offset.epoch(), target);
        assert_state(
            &at,
            offset.position().as_array().map(|x| x.value()),
            offset.velocity().as_array().map(|x| x.value()),
            1e-9,
            1e-12,
        );
    }
    let leap = Time::<UTC>::try_from_chrono_with("2016-12-31T23:59:60Z".parse().unwrap(), &ctx)
        .unwrap()
        .reinterpret::<JD>();
    assert!(matches!(
        prop.propagate_at(leap),
        Err(Sgp4Error::TimeConversion(_))
    ));
    tle.epoch = leap.reinterpret();
    assert!(Sgp4Propagator::from_tle(&tle).is_err());
}

#[test]
fn iau_models_match_independent_vallado_cpp_port() {
    let mut lines = include_str!("fixtures/sgp4/iau.txt")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .peekable();
    let mut checked = 0;
    while let Some(line1) = lines.next() {
        let tle = parse_tle(line1, lines.next().unwrap()).unwrap();
        while lines.peek().is_some_and(|l| !l.starts_with("1 ")) {
            let values: Vec<_> = lines.next().unwrap().split_whitespace().collect();
            let model = match values[0] {
                "72" => GravityModel::Wgs72Iau,
                "84" => GravityModel::Wgs84,
                _ => panic!("unknown fixture model"),
            };
            let prop = Sgp4Propagator::from_tle_with_model(&tle, model).unwrap();
            assert_eq!(prop.gravity_model(), model);
            let numbers: Vec<f64> = values[1..].iter().map(|s| s.parse().unwrap()).collect();
            let state = prop.propagate_minutes(Minutes::new(numbers[0])).unwrap();
            assert_state(
                &state,
                numbers[1..4].try_into().unwrap(),
                numbers[4..].try_into().unwrap(),
                1e-6,
                1e-9,
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 48);
}

#[test]
fn metadata_does_not_change_propagation_and_requests_are_order_independent() {
    use siderust::formats::tle::Classification;
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Sgp4Propagator>();
    let valid = parse_tle(L1, L2).unwrap();
    let prop = Sgp4Propagator::from_tle(&valid).unwrap();
    let expected = prop.propagate_minutes(Minutes::new(-1440.0)).unwrap();
    prop.propagate_minutes(Minutes::new(1440.0)).unwrap();
    for class in [
        Classification::Unclassified,
        Classification::Classified,
        Classification::Secret,
    ] {
        let mut tle = valid.clone();
        tle.classification = class;
        tle.mean_motion_dot = 1.0;
        tle.mean_motion_ddot = -2.0;
        let copy = Sgp4Propagator::from_tle(&tle).unwrap().clone();
        let state = copy.propagate_minutes(Minutes::new(-1440.0)).unwrap();
        assert_eq!(state.position(), expected.position());
        assert_eq!(state.velocity(), expected.velocity());
    }
}

#[test]
fn extreme_finite_elements_cannot_return_non_finite_states() {
    let mut tle = parse_tle(L1, L2).unwrap();
    tle.bstar = f64::MAX;
    assert!(matches!(
        Sgp4Propagator::from_tle(&tle),
        Err(Sgp4Error::InvalidElements { .. })
    ));
    tle.bstar = 1e5;
    let prop = Sgp4Propagator::from_tle(&tle).unwrap();
    assert!(matches!(
        prop.propagate_minutes(Minutes::new(10_000.0)),
        Err(Sgp4Error::Propagation { .. })
    ));

    // Circular low-perigee orbit: the eccentricity drag vanishes, but the
    // deliberately extreme finite BSTAR overflows the radius calculation.
    // This must reach the Cartesian-output guard rather than return NaNs.
    tle.bstar = f64::MAX;
    tle.eccentricity = 0.0;
    tle.inclination = Degrees::new(54.735610317245346);
    tle.mean_motion = qtty::Quantity::new(16.5);
    let prop = Sgp4Propagator::from_tle(&tle).unwrap();
    let error = prop.propagate_minutes(Minutes::new(1.0)).unwrap_err();
    assert!(
        error.to_string().contains("non-finite Cartesian state"),
        "{error}"
    );
}

#[test]
fn fractional_minutes_preserve_precision_through_typed_epochs() {
    let prop = Sgp4Propagator::from_tle(&parse_tle(L1, L2).unwrap()).unwrap();
    for minutes in [
        -1440.123456789,
        -0.00000001,
        0.00000001,
        0.123456789,
        1440.123456789,
    ] {
        let offset = prop.propagate_minutes(Minutes::new(minutes)).unwrap();
        let at = prop.propagate_at(offset.epoch()).unwrap();
        assert_state(
            &at,
            offset.position().as_array().map(|q| q.value()),
            offset.velocity().as_array().map(|q| q.value()),
            1e-6,
            1e-9,
        );
    }
    let distant = Time::<UTC>::try_from_chrono_with(
        "2400-01-01T00:00:00Z".parse().unwrap(),
        prop.time_context(),
    )
    .unwrap()
    .to_with::<JD>(prop.time_context())
    .unwrap();
    assert!(matches!(
        prop.propagate_at(distant),
        Err(Sgp4Error::TimeConversion(_))
    ));
    assert!(matches!(
        prop.propagate_minutes(Minutes::new(-100.0 * 365.25 * 1440.0)),
        Err(Sgp4Error::TimeConversion(_))
    ));
}

#[test]
fn synthetic_boundary_cases_match_independent_cpp() {
    let mut tle = parse_tle(L1, L2).unwrap();
    let ctx = TimeContext::new();
    tle.epoch =
        Time::<UTC>::try_from_chrono_with("2020-01-01T00:00:00Z".parse().unwrap(), &ctx).unwrap();
    tle.raan = Degrees::new(-14.0);
    tle.argument_of_perigee = Degrees::new(41.0);
    tle.mean_anomaly = Degrees::new(17.0);
    let mut checked = 0;
    for line in include_str!("fixtures/sgp4/synthetic.txt")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split_whitespace().collect();
        let model = match fields[0] {
            "72a" => GravityModel::Wgs72,
            "72i" => GravityModel::Wgs72Iau,
            "84i" => GravityModel::Wgs84,
            _ => panic!("unknown model"),
        };
        let v: Vec<f64> = fields[1..].iter().map(|x| x.parse().unwrap()).collect();
        tle.eccentricity = v[0];
        tle.inclination = Degrees::new(v[1]);
        tle.mean_motion = qtty::Quantity::new(v[2]);
        tle.bstar = v[3];
        let prop = Sgp4Propagator::from_tle_with_model(&tle, model).unwrap();
        let state = prop.propagate_minutes(Minutes::new(v[4])).unwrap();
        assert_state(
            &state,
            v[5..8].try_into().unwrap(),
            v[8..].try_into().unwrap(),
            1e-6,
            1e-9,
        );
        checked += 1;
    }
    assert_eq!(checked, 315);
}
