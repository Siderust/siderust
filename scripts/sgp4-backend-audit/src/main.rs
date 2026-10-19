// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Reproduce upstream parity failures against unchanged independent C++ vectors.
//! This executable deliberately exits unsuccessfully when parity fails.
use std::f64::consts::PI;

fn main() {
    let mut counts = [[0; 2]; 3];
    let mut maxima = [[0.0_f64; 2]; 3];
    for (index, line) in include_str!("../../../tests/fixtures/sgp4/synthetic.txt")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .enumerate()
    {
        let fields: Vec<_> = line.split_whitespace().collect();
        let model = match fields[0] {
            "72a" => 0,
            "72i" => 1,
            "84i" => 2,
            _ => panic!("unknown fixture gravity model"),
        };
        let v: Vec<f64> = fields[1..].iter().map(|s| s.parse().unwrap()).collect();
        let geopotential = if model == 2 { sgp4::WGS84 } else { sgp4::WGS72 };
        let orbit = sgp4::Orbit::from_kozai_elements(
            &geopotential,
            v[1].to_radians(),
            (-14.0_f64).to_radians(),
            v[0],
            41.0_f64.to_radians(),
            17.0_f64.to_radians(),
            v[2] * PI / 720.0,
        )
        .unwrap();
        // Fixture epoch: 2020-01-01 00:00 UTC, exactly 7304.5 days after J2000.
        // An integer half-day is represented exactly in both epoch conventions.
        let constants = sgp4::Constants::new(
            geopotential,
            if model == 0 {
                sgp4::afspc_epoch_to_sidereal_time
            } else {
                sgp4::iau_epoch_to_sidereal_time
            },
            7304.5 / 365.25,
            v[3],
            orbit,
        )
        .unwrap();
        let offset = sgp4::MinutesSinceEpoch(v[4]);
        let prediction = if model == 0 {
            constants.propagate_afspc_compatibility_mode(offset)
        } else {
            constants.propagate(offset)
        }
        .unwrap();
        let position_error = prediction
            .position
            .iter()
            .zip(&v[5..8])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let velocity_error = prediction
            .velocity
            .iter()
            .zip(&v[8..])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        counts[model][0] += 1;
        maxima[model][0] = maxima[model][0].max(position_error);
        maxima[model][1] = maxima[model][1].max(velocity_error);
        let non_finite = prediction
            .position
            .iter()
            .chain(&prediction.velocity)
            .any(|value| !value.is_finite());
        if non_finite
            || !position_error.is_finite()
            || !velocity_error.is_finite()
            || position_error >= 1e-6
            || velocity_error >= 1e-9
        {
            counts[model][1] += 1;
            println!("FAIL row={} model={} e={} i={} minutes={} position={position_error:e} km velocity={velocity_error:e} km/s",
                index + 1, fields[0], v[0], v[1], v[4]);
        }
    }
    for (i, name) in ["Wgs72", "Wgs72Iau", "Wgs84"].iter().enumerate() {
        println!(
            "{name}: {}/{} failed; max position={} km; max velocity={} km/s",
            counts[i][1], counts[i][0], maxima[i][0], maxima[i][1]
        );
    }
    if counts.iter().any(|count| count[1] > 0) {
        std::process::exit(1);
    }
}
