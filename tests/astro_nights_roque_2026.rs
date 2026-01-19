use std::fs::File;
use std::io::BufReader;

use chrono::{NaiveDate, NaiveDateTime, NaiveTime, TimeZone, Utc};
use siderust::calculus::solar::altitude_periods::{find_night_periods, twilight};
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};

/// Compares siderust-computed astronomical nights for Roque de los Muchachos
/// against the JSON file placed at `tests/reference_data/astro_night_periods_2026.json`.
#[test]
fn test_astronomical_nights_roque_2026() {
    let f = File::open("tests/reference_data/astro_night_periods_2026.json")
        .expect("Missing tests/reference_data/astro_night_periods_2026.json");
    let reader = BufReader::new(f);
    let json: serde_json::Value =
        serde_json::from_reader(reader).expect("Invalid JSON in reference file");
    let expected_periods: Vec<Period<ModifiedJulianDate>> =
        serde_json::from_value(json["periods"].clone())
            .expect("Invalid or missing `periods` in reference file");

    // Build observer site from catalog constant
    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    // Build search period: 2026-01-01T00:00:00Z -> 2027-01-01T00:00:00Z
    let start_naive = NaiveDate::from_ymd_opt(2026, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(2027, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());

    let start_dt =
        Utc.from_utc_datetime(&NaiveDateTime::new(start_naive.date(), start_naive.time()));
    let end_dt = Utc.from_utc_datetime(&NaiveDateTime::new(end_naive.date(), end_naive.time()));

    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);
    let period = Period::new(mjd_start, mjd_end);

    let computed_opt = find_night_periods(site, period, twilight::ASTRONOMICAL);
    assert!(
        computed_opt.is_some(),
        "siderust did not find any astronomical nights"
    );
    let computed = computed_opt.unwrap();

    assert_eq!(
        expected_periods.len(),
        computed.len(),
        "Different number of intervals between expected and siderust computation"
    );

    // Tolerance: 3 seconds (expressed in days)
    let tol_seconds = 3.0_f64;
    let tol_days = tol_seconds / 86400.0;

    for (i, (exp_p, calc_p)) in expected_periods.iter().zip(computed.iter()).enumerate() {
        let ds = (exp_p.start.value() - calc_p.start.value()).abs();
        let de = (exp_p.end.value() - calc_p.end.value()).abs();

        assert!(
            ds <= tol_days,
            "Start MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.start.value(),
            calc_p.start.value(),
            ds,
            ds * 86400.0
        );

        assert!(
            de <= tol_days,
            "End MJD mismatch at index {}: expected {} got {} (diff {} days ~ {} s)",
            i,
            exp_p.end.value(),
            calc_p.end.value(),
            de,
            de * 86400.0
        );
    }
}
