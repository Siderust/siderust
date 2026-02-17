// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![cfg(feature = "serde")]

use std::fs::File;
use std::io::BufReader;

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use serde_json::Value;
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period, UniversalTime, MJD};

const REFERENCE_PATH: &str = "tests/reference_data/astro_night_periods_2026.json";
const TOL_SECONDS: f64 = 3.0;
const TOL_DAYS: f64 = TOL_SECONDS / 86_400.0;

fn load_reference_periods() -> Vec<Period<MJD>> {
    let f = File::open(REFERENCE_PATH).expect("Missing astro_night_periods_2026.json");
    let reader = BufReader::new(f);
    let json: Value = serde_json::from_reader(reader).expect("Invalid JSON in reference file");
    let periods_utc_axis: Vec<Period<MJD>> = serde_json::from_value(json["periods"].clone())
        .expect("Invalid or missing `periods` in reference file");

    // JSON stores UTC-based MJD values; this API expects MJD values on the TT axis.
    periods_utc_axis
        .into_iter()
        .map(|p| {
            Period::new(
                mjd_utc_axis_to_mjd_tt(p.start.value()),
                mjd_utc_axis_to_mjd_tt(p.end.value()),
            )
        })
        .collect()
}

fn mjd_utc_axis_to_mjd_tt(mjd_utc: f64) -> ModifiedJulianDate {
    // tempoch uses JD(TT) as canonical axis. Convert JD(UT/UTC-like) -> JD(TT).
    const MJD_OFFSET: f64 = 2_400_000.5;
    UniversalTime::new(mjd_utc + MJD_OFFSET).to::<MJD>()
}

fn build_roque_period() -> (ObserverSite, Period<MJD>) {
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(2027, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    // Query window endpoints are expressed in UTC, then converted to the TT axis.
    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);

    (site, Period::new(mjd_start, mjd_end))
}

fn assert_periods_close(expected: &[Period<MJD>], computed: &[Period<MJD>]) {
    assert_eq!(
        expected.len(),
        computed.len(),
        "Different number of intervals between expected and siderust computation"
    );

    for (i, (exp_p, calc_p)) in expected.iter().zip(computed.iter()).enumerate() {
        let ds = (exp_p.start.value() - calc_p.start.value()).abs();
        let de = (exp_p.end.value() - calc_p.end.value()).abs();

        assert!(
            ds <= TOL_DAYS,
            "Start MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.start.value(),
            calc_p.start.value(),
            ds,
            ds * 86_400.0
        );

        assert!(
            de <= TOL_DAYS,
            "End MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.end.value(),
            calc_p.end.value(),
            de,
            de * 86_400.0
        );
    }
}

/// Compares siderust-computed astronomical nights for Roque de los Muchachos
/// against the JSON reference data.
#[test]
fn test_astronomical_nights_roque_2026() {
    let expected_periods = load_reference_periods();
    let (site, period) = build_roque_period();

    let computed = Sun.below_threshold(site, period, twilight::ASTRONOMICAL);
    assert!(
        !computed.is_empty(),
        "siderust did not find any astronomical nights"
    );

    assert_periods_close(&expected_periods, &computed);
}

/// Same as the default test but exercises the culmination-based altitude search
/// followed by complement to obtain night periods.
#[test]
fn test_astronomical_nights_roque_2026_culminations() {
    let expected_periods = load_reference_periods();
    let (site, period) = build_roque_period();

    // above_threshold finds "above" periods; complement gives us the nights.
    let day_periods = Sun.above_threshold(site, period, twilight::ASTRONOMICAL);
    let computed = siderust::time::complement_within(period, &day_periods);
    assert!(
        !computed.is_empty(),
        "Culmination-based search did not find any astronomical nights"
    );

    assert_periods_close(&expected_periods, &computed);
}
