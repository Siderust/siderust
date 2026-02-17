// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![cfg(feature = "serde")]

use std::fs::File;
use std::io::BufReader;

use chrono::{DateTime, Utc};
use serde::Deserialize;
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period, MJD};

const REFERENCE_PATH: &str = "tests/reference_data/astro_night_periods_2026.json";
const TOL_SECONDS: f64 = 3.0;
const TOL_DAYS: f64 = TOL_SECONDS / 86_400.0;
const UNIX_EPOCH_MJD_UTC: f64 = 40_587.0;
const SECONDS_PER_DAY: f64 = 86_400.0;

#[derive(Debug, Deserialize)]
struct ReferencePeriodUtcMjd {
    start_mjd: f64,
    end_mjd: f64,
}

#[derive(Debug, Deserialize)]
struct ReferenceDataFile {
    start_utc: String,
    end_utc: String,
    period_mjd_scale: String,
    periods: Vec<ReferencePeriodUtcMjd>,
}

#[derive(Debug)]
struct ReferenceData {
    window_tt: Period<MJD>,
    periods_tt: Vec<Period<MJD>>,
}

fn parse_utc_datetime(s: &str) -> DateTime<Utc> {
    DateTime::parse_from_rfc3339(s)
        .unwrap_or_else(|err| panic!("Invalid RFC3339 timestamp `{s}`: {err}"))
        .with_timezone(&Utc)
}

fn utc_mjd_to_utc_datetime(mjd_utc: f64) -> DateTime<Utc> {
    let total_seconds = (mjd_utc - UNIX_EPOCH_MJD_UTC) * SECONDS_PER_DAY;
    let mut whole_seconds = total_seconds.floor() as i64;
    let mut nanos = ((total_seconds - whole_seconds as f64) * 1_000_000_000.0).round() as i64;

    if nanos == 1_000_000_000 {
        whole_seconds += 1;
        nanos = 0;
    } else if nanos < 0 {
        whole_seconds -= 1;
        nanos += 1_000_000_000;
    }

    DateTime::<Utc>::from_timestamp(whole_seconds, nanos as u32)
        .unwrap_or_else(|| panic!("UTC MJD `{mjd_utc}` is out of DateTime<Utc> range"))
}

fn utc_mjd_to_tt_mjd(mjd_utc: f64) -> ModifiedJulianDate {
    ModifiedJulianDate::from_utc(utc_mjd_to_utc_datetime(mjd_utc))
}

fn load_reference_data() -> ReferenceData {
    let f = File::open(REFERENCE_PATH).expect("Missing astro_night_periods_2026.json");
    let reader = BufReader::new(f);
    let reference: ReferenceDataFile =
        serde_json::from_reader(reader).expect("Invalid JSON in reference file");

    assert_eq!(
        reference.period_mjd_scale.as_str(),
        "UTC",
        "Reference data period MJD values must be tagged as UTC scale"
    );

    let start_utc = parse_utc_datetime(&reference.start_utc);
    let end_utc = parse_utc_datetime(&reference.end_utc);
    let window_tt = Period::new(
        ModifiedJulianDate::from_utc(start_utc),
        ModifiedJulianDate::from_utc(end_utc),
    );

    // Reference JSON stores period endpoints as UTC-based MJD numbers.
    // Siderust altitude APIs use Period<MJD> values on the canonical TT axis.
    let periods_tt = reference
        .periods
        .into_iter()
        .map(|period| {
            Period::new(
                utc_mjd_to_tt_mjd(period.start_mjd),
                utc_mjd_to_tt_mjd(period.end_mjd),
            )
        })
        .collect();

    ReferenceData {
        window_tt,
        periods_tt,
    }
}

fn roque_site() -> ObserverSite {
    ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS)
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
            "Start TT-MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.start.value(),
            calc_p.start.value(),
            ds,
            ds * 86_400.0
        );

        assert!(
            de <= TOL_DAYS,
            "End TT-MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.end.value(),
            calc_p.end.value(),
            de,
            de * 86_400.0
        );
    }
}

/// Compares siderust-computed astronomical nights for Roque de los Muchachos
/// against JSON reference data whose period MJD values are tagged as UTC and
/// converted to TT MJD for comparison.
#[test]
fn test_astronomical_nights_roque_2026() {
    let reference = load_reference_data();
    let site = roque_site();

    let computed = Sun.below_threshold(site, reference.window_tt, twilight::ASTRONOMICAL);
    assert!(
        !computed.is_empty(),
        "siderust did not find any astronomical nights"
    );

    assert_periods_close(&reference.periods_tt, &computed);
}

/// Same as the default test but exercises the culmination-based altitude search
/// followed by complement to obtain night periods.
#[test]
fn test_astronomical_nights_roque_2026_culminations() {
    let reference = load_reference_data();
    let site = roque_site();

    // above_threshold finds "above" periods; complement gives us the nights.
    let day_periods = Sun.above_threshold(site, reference.window_tt, twilight::ASTRONOMICAL);
    let computed = siderust::time::complement_within(reference.window_tt, &day_periods);
    assert!(
        !computed.is_empty(),
        "Culmination-based search did not find any astronomical nights"
    );

    assert_periods_close(&reference.periods_tt, &computed);
}
