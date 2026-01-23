// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Benchmarks for lunar altitude calculations.
//!
//! Tests the performance of calculating Moon altitude and finding periods when
//! the Moon is above/below the horizon using the ELP2000 lunar theory.

use criterion::{criterion_group, criterion_main, Criterion};
use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use qtty::*;
use siderust::astro::JulianDate;
use siderust::bodies::solar_system::Moon;
use siderust::calculus::events::altitude_periods::{find_altitude_periods, AltitudeCondition};
use siderust::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use siderust::coordinates::frames::{self, Ecliptic};
use siderust::coordinates::spherical::Position;
use siderust::coordinates::transform::{Transform, TransformFrame};
use siderust::coordinates::transform::centers::ToTopocentricExt;
use siderust::coordinates::cartesian;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};
use std::hint::black_box;
use std::time::Duration;

/// Compute Moon altitude at a given Julian Date for a specific observer site.
/// This is the core function being benchmarked.
fn moon_altitude_rad(jd: JulianDate, site: &ObserverSite) -> f64 {
    // Get Moon's geocentric ecliptic position from ELP2000 (cartesian)
    let moon_geo_ecliptic_cart: cartesian::Position<Geocentric, Ecliptic, Kilometer> =
        Moon::get_geo_position(jd);

    // Transform: Ecliptic -> EquatorialMeanJ2000 (cartesian)
    let moon_geo_eq_j2000_cart: cartesian::Position<
        Geocentric,
        frames::EquatorialMeanJ2000,
        Kilometer,
    > = TransformFrame::to_frame(&moon_geo_ecliptic_cart);

    // Apply topocentric parallax correction (important for the Moon!)
    let moon_topo_eq_j2000_cart: cartesian::Position<
        Topocentric,
        frames::EquatorialMeanJ2000,
        Kilometer,
    > = moon_geo_eq_j2000_cart.to_topocentric(*site, jd);

    // Apply precession: J2000 -> mean-of-date
    let rot = siderust::coordinates::transform::frame_rotation::<
        frames::EquatorialMeanJ2000,
        frames::EquatorialMeanOfDate,
    >(jd, &siderust::coordinates::transform::AstroContext::default());

    let [x, y, z] = rot.apply_array([
        moon_topo_eq_j2000_cart.x().value(),
        moon_topo_eq_j2000_cart.y().value(),
        moon_topo_eq_j2000_cart.z().value(),
    ]);

    let moon_topo_eq_date: cartesian::Position<
        Topocentric,
        frames::EquatorialMeanOfDate,
        Kilometer,
    > = cartesian::Position::from_vec3(
        *site,
        nalgebra::Vector3::new(x * KM, y * KM, z * KM),
    );

    // Transform to horizontal coordinates
    let moon_horizontal: cartesian::Position<Topocentric, frames::Horizontal, Kilometer> =
        moon_topo_eq_date.transform(jd);

    // Convert to spherical to extract altitude
    let moon_horizontal_spherical: Position<Topocentric, frames::Horizontal, Kilometer> =
        Position::<Topocentric, frames::Horizontal, Kilometer>::from_cartesian(&moon_horizontal);

    moon_horizontal_spherical.alt().to::<Radian>().value()
}

fn build_period(days: u32) -> Period<ModifiedJulianDate> {
    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap())
        + chrono::Duration::days(days as i64);

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);

    Period::new(mjd_start, mjd_end)
}

fn bench_moon_altitude_computation(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    let jd = JulianDate::J2000;

    let mut group = c.benchmark_group("moon_altitude_single");

    // Benchmark single altitude computation
    group.bench_function("compute_altitude", |b| {
        b.iter(|| {
            let _altitude = moon_altitude_rad(black_box(jd), black_box(&site));
        });
    });

    group.finish();
}

fn bench_moon_above_horizon(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_above_horizon");

    // Benchmark for 1-day horizon
    group.bench_function("find_moon_above_horizon_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::above(Degrees::new(0.0))),
            );
        });
    });

    // Benchmark for 7-day horizon
    group.bench_function("find_moon_above_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::above(Degrees::new(0.0))),
            );
        });
    });

    // Benchmark for 30-day horizon (full lunar cycle)
    group.bench_function("find_moon_above_horizon_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::above(Degrees::new(0.0))),
            );
        });
    });

    group.finish();
}

fn bench_moon_below_horizon(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_below_horizon");

    // Benchmark for 1-day horizon
    group.bench_function("find_moon_below_horizon_1day", |b| {
        let period = black_box(build_period(1));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::below(Degrees::new(-0.5))),
            );
        });
    });

    // Benchmark for 7-day horizon
    group.bench_function("find_moon_below_horizon_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::below(Degrees::new(-0.5))),
            );
        });
    });

    // Benchmark for 30-day horizon
    group.bench_function("find_moon_below_horizon_30day", |b| {
        let period = black_box(build_period(30));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::below(Degrees::new(-0.5))),
            );
        });
    });

    group.finish();
}

fn bench_moon_altitude_range(c: &mut Criterion) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let mut group = c.benchmark_group("moon_altitude_range");

    // Benchmark for finding Moon at low altitude (0-30 degrees) over 7 days
    group.bench_function("find_moon_low_altitude_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::Between {
                    min: Degrees::new(0.0),
                    max: Degrees::new(30.0),
                }),
            );
        });
    });

    // Benchmark for finding Moon at high altitude (60-90 degrees) over 7 days
    group.bench_function("find_moon_high_altitude_7day", |b| {
        let period = black_box(build_period(7));
        b.iter(|| {
            let altitude_fn = |jd: JulianDate| -> f64 { moon_altitude_rad(jd, &site) };
            let _result = find_altitude_periods(
                black_box(altitude_fn),
                black_box(period),
                black_box(AltitudeCondition::Between {
                    min: Degrees::new(60.0),
                    max: Degrees::new(90.0),
                }),
            );
        });
    });

    group.finish();
}

criterion_group! {
    name = moon_benches;
    config = Criterion::default()
        .measurement_time(Duration::from_secs(10))
        .sample_size(20);
    targets = bench_moon_altitude_computation, bench_moon_above_horizon, bench_moon_below_horizon, bench_moon_altitude_range
}
criterion_main!(moon_benches);
