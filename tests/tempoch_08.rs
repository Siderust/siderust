//! Regression values for the civil UTC integration introduced with tempoch 0.8.
use chrono::{DateTime, Timelike, Utc};
use qtty::time::Minutes;
use qtty::Day;
use siderust::astro::sgp4::TlePropagator;
use siderust::formats::tle::parse_tle;
use tempoch::{
    ConversionError, JulianDate, ModifiedJulianDate, Time, TimeContext, JD, MJD, TAI, TT, UTC,
};

fn calendar(label: &str, ctx: &TimeContext) -> Time<UTC> {
    let dt: DateTime<Utc> = label.parse().unwrap();
    Time::<UTC>::try_from_chrono_with(dt, ctx).unwrap()
}

fn assert_calendar(actual: DateTime<Utc>, label: &str) {
    let expected: DateTime<Utc> = label.parse().unwrap();
    let error = actual
        .signed_duration_since(expected)
        .num_nanoseconds()
        .unwrap();
    assert!(error.abs() < 100_000, "civil day scalar error: {error} ns");
}

#[test]
fn civil_mjd_2017_has_known_calendar_and_atomic_offsets() {
    let ctx = TimeContext::new();
    let utc = ModifiedJulianDate::<UTC>::try_new_with(Day::new(57754.0), &ctx).unwrap();
    assert_calendar(
        utc.try_to_chrono_with(&ctx).unwrap(),
        "2017-01-01T00:00:00+00:00",
    );
    assert_eq!(
        utc.to_with::<JD>(&ctx)
            .unwrap()
            .try_raw_with(&ctx)
            .unwrap()
            .value(),
        2457754.5
    );
    let tai = utc.to::<TAI>().to::<MJD>().raw().value();
    let tt = utc.to::<TT>().to::<MJD>().raw().value();
    assert!(((tai - 57754.0) * 86400.0 - 37.0).abs() < 1e-4);
    assert!(((tt - 57754.0) * 86400.0 - 69.184).abs() < 1e-4);
    // Continuous formats retain infallible construction and extraction.
    assert_eq!(
        JulianDate::<TT>::new(2451545.0).to::<MJD>().raw().value(),
        51544.5
    );
    assert_eq!(
        JulianDate::<TAI>::new(2451545.0).to::<MJD>().raw().value(),
        51544.5
    );
}

#[test]
fn civil_days_and_si_days_differ_at_the_2016_leap_second() {
    let ctx = TimeContext::new();
    let before = ModifiedJulianDate::<UTC>::try_new_with(Day::new(57753.0), &ctx).unwrap();
    let after = ModifiedJulianDate::<UTC>::try_new_with(Day::new(57754.0), &ctx).unwrap();
    assert!((after.to_j2000s().duration_since(before.to_j2000s()).value() - 86401.0).abs() < 1e-9);
    let leap = calendar("2016-12-31T23:59:60Z", &ctx);
    assert_eq!(
        leap.to_with::<JD>(&ctx).unwrap_err(),
        ConversionError::InvalidLeapSecond
    );
    assert_eq!(
        leap.to_with::<MJD>(&ctx).unwrap_err(),
        ConversionError::InvalidLeapSecond
    );
    // Reinterpreting a phantom format cannot bypass upstream validation.
    assert_eq!(
        leap.reinterpret::<MJD>().try_raw_with(&ctx).unwrap_err(),
        ConversionError::InvalidLeapSecond
    );
}

#[test]
fn tle_model_uses_nominal_minutes_across_leap_midnight() {
    let mut tle = parse_tle(
        "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
        "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
    )
    .unwrap();
    let ctx = TimeContext::new();
    tle.epoch = calendar("2016-12-31T00:00:00Z", &ctx);
    let prop = TlePropagator::from_tle(&tle).unwrap();
    assert_eq!(
        prop.epoch_jd_utc()
            .try_raw_with(prop.time_context())
            .unwrap()
            .value(),
        2457753.5
    );
    let next_midnight = calendar("2017-01-01T00:00:00Z", &ctx)
        .to_with::<JD>(&ctx)
        .unwrap();
    let at = prop.propagate_at(next_midnight).unwrap();
    let minutes = prop.propagate_minutes(Minutes::new(1440.0)).unwrap();
    assert_eq!(at.position(), minutes.position());
    assert_eq!(at.velocity(), minutes.velocity());
    assert_calendar(
        minutes.epoch().try_to_chrono_with(&ctx).unwrap(),
        "2017-01-01T00:00:00Z",
    );
    assert!(
        (minutes
            .epoch()
            .to_j2000s()
            .duration_since(prop.epoch_jd_utc().to_j2000s())
            .value()
            - 86401.0)
            .abs()
            < 1e-4
    );
    let leap = calendar("2016-12-31T23:59:60Z", &ctx).reinterpret::<JD>();
    assert!(prop.propagate_at(leap).is_err());
    assert!(prop.propagate_minutes(Minutes::new(f64::MAX)).is_err());
}

#[cfg(feature = "std")]
#[test]
fn c04_interpolation_follows_civil_mjd_not_elapsed_si_time() {
    use siderust::formats::iers::eop::{interpolate, read_eop_c04};
    let records = read_eop_c04(
        "2016 12 31 57753 0.1 0.2 -0.4 0.0 0.0 0.0\n2017 1 1 57754 0.3 0.4 0.6 0.0 0.0 0.0\n"
            .as_bytes(),
    )
    .unwrap();
    let ctx = TimeContext::new();
    assert_calendar(
        records[0].mjd.try_to_chrono_with(&ctx).unwrap(),
        "2016-12-31T00:00:00+00:00",
    );
    let noon = calendar("2016-12-31T12:00:00Z", &ctx)
        .to_with::<MJD>(&ctx)
        .unwrap();
    let mid = interpolate(&records, noon, &ctx).unwrap().unwrap();
    assert!((mid.eop.xp.value() - 0.2).abs() < 1e-8);
    assert!((mid.eop.dut1.value() + 0.4).abs() < 1e-8);
    let leap = calendar("2016-12-31T23:59:60Z", &ctx).reinterpret::<MJD>();
    assert_eq!(
        interpolate(&records, leap, &ctx).unwrap_err(),
        ConversionError::InvalidLeapSecond
    );
}

#[cfg(feature = "std")]
#[test]
fn cpf_epoch_is_civil_mjd_with_fractional_seconds() {
    use siderust::formats::ilrs::cpf::parse_cpf;
    let cpf = parse_cpf("H1 CPF 2 TEST 2017 01 01 00 1\nH2 TEST 1 5 2017 01 01 00 00 00 2017 01 02 00 00 00 60 0 0\n10 1 57754 43200.125 0 1000000 2000000 3000000\n99\n").unwrap();
    let epoch = cpf.ephemeris.entries[0].epoch.try_to_chrono().unwrap();
    let expected: DateTime<Utc> = "2017-01-01T12:00:00.125Z".parse().unwrap();
    assert!(
        (epoch
            .signed_duration_since(expected)
            .num_nanoseconds()
            .unwrap())
        .abs()
            < 100_000
    );
}

#[cfg(feature = "serde")]
#[test]
fn serde_uses_civil_utc_scalars_and_rejects_leap_instants() {
    let ctx = TimeContext::new();
    let epoch: ModifiedJulianDate<UTC> = serde_json::from_str("57754.0").unwrap();
    assert_calendar(
        epoch.try_to_chrono_with(&ctx).unwrap(),
        "2017-01-01T00:00:00+00:00",
    );
    assert!(
        (serde_json::to_string(&epoch)
            .unwrap()
            .parse::<f64>()
            .unwrap()
            - 57754.0)
            .abs()
            < 1e-9
    );
    let leap = calendar("2016-12-31T23:59:60Z", &ctx).reinterpret::<MJD>();
    assert!(serde_json::to_string(&leap).is_err());
}

#[test]
fn utc_context_policy_and_invalid_coordinates_propagate_errors() {
    let ctx = TimeContext::new();
    assert_eq!(
        ModifiedJulianDate::<UTC>::try_new_with(Day::new(0.0), &ctx).unwrap_err(),
        ConversionError::UtcBeforeDefinition
    );
    let historical = ctx.clone().allow_pre_definition_utc();
    assert!(ModifiedJulianDate::<UTC>::try_new_with(Day::new(0.0), &historical).is_ok());
    for raw in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            ModifiedJulianDate::<UTC>::try_new_with(Day::new(raw), &ctx).unwrap_err(),
            ConversionError::NonFinite
        );
        assert_eq!(
            JulianDate::<UTC>::try_new_with(Day::new(raw), &ctx).unwrap_err(),
            ConversionError::NonFinite
        );
    }
}

#[cfg(feature = "std")]
#[test]
fn cpf_preserves_leap_second_labels_and_fractional_seconds() {
    use siderust::formats::ilrs::cpf::{parse_cpf, write_cpf};
    let ctx = TimeContext::new();
    for sod in [86400.0, 86400.125, 86400.999, 86400.999999999] {
        let input = format!("H1 CPF 2 TEST 2016 12 31 00 1\nH2 TEST 1 5 2016 12 31 00 00 00 2017 01 01 00 00 00 60 0 0\n10 1 57753 {sod:.9} 0 1000000 2000000 3000000\n99\n");
        let cpf = parse_cpf(&input).unwrap();
        let civil = cpf.ephemeris.entries[0]
            .epoch
            .try_to_chrono_with(&ctx)
            .unwrap();
        assert_eq!(civil.second(), 59);
        assert!(
            (cpf.ephemeris.entries[0]
                .epoch
                .duration_since(
                    ModifiedJulianDate::<UTC>::try_new_with(Day::new(57753.0), &ctx)
                        .unwrap()
                        .to_j2000s()
                )
                .value()
                - sod)
                .abs()
                < 2e-7
        );
        let mut output = Vec::new();
        write_cpf(&cpf, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains(&format!("10 1 57753 {sod:.9}")), "{output}");
        let roundtrip = parse_cpf(&output).unwrap();
        assert!(
            (roundtrip.ephemeris.entries[0]
                .epoch
                .duration_since(cpf.ephemeris.entries[0].epoch)
                .value())
            .abs()
                < 2e-7
        );
        assert!(parse_cpf(&input.replace("10 1 57753", "10 1 57754")).is_err());
    }
}

#[cfg(feature = "std")]
#[test]
fn eop_dut1_jump_is_preserved_at_leap_boundary() {
    use siderust::formats::iers::eop::{interpolate, read_eop_c04};
    let records = read_eop_c04("2016 12 31 57753 0.1 0.2 -0.4 0.0 0.0 0.0\n2017 1 1 57754 0.3 0.4 0.6 0.0 0.0 0.0\n2017 1 2 57755 0.5 0.6 0.4 0.0 0.0 0.0\n".as_bytes()).unwrap();
    let ctx = TimeContext::new();
    for (label, dut1) in [
        ("2016-12-31T23:59:59Z", -0.4),
        ("2017-01-01T00:00:00Z", 0.6),
        ("2017-01-01T12:00:00Z", 0.5),
    ] {
        let query = calendar(label, &ctx).to_with::<MJD>(&ctx).unwrap();
        let value = interpolate(&records, query, &ctx).unwrap().unwrap();
        assert!(
            (value.eop.dut1.value() - dut1).abs() < 2e-7,
            "{label}: {}",
            value.eop.dut1
        );
    }
}

#[test]
fn earth_rotation_helpers_preserve_ut1_scale_and_supplied_dut1() {
    use siderust::astro::earth_rotation::{jd_ut1_from_tt, jd_ut1_from_tt_eop};
    let ctx = TimeContext::new();
    let utc = calendar("2016-12-31T23:59:60.125Z", &ctx);
    let tt = utc.to::<TT>().to::<JD>();
    let _: JulianDate<tempoch::UT1> = jd_ut1_from_tt(tt);
    let eop = siderust::astro::eop::EopValues {
        dut1: qtty::time::Seconds::new(-0.4),
        ..Default::default()
    };
    let ut1: JulianDate<tempoch::UT1> = jd_ut1_from_tt_eop(tt, &eop);
    let expected = utc.to::<TAI>().to_j2000s() - qtty::Second::new(36.4);
    let (actual_hi, actual_lo) = ut1.raw_seconds_pair();
    let (expected_hi, expected_lo) = expected.raw_seconds_pair();
    assert!(
        ((actual_hi - expected_hi) + (actual_lo - expected_lo))
            .value()
            .abs()
            < 1e-8
    );
}

#[cfg(feature = "std")]
#[test]
fn invalid_c04_values_report_errors() {
    use siderust::formats::iers::eop::read_eop_c04;
    for value in ["invalid", "NaN", "inf"] {
        let line = format!("2017 1 1 57754 0.1 0.2 {value} 0.0 0.0 0.0\n");
        assert!(read_eop_c04(line.as_bytes()).is_err());
    }
}

#[cfg(feature = "std")]
#[test]
fn cpf_roundtrips_the_civil_labels_adjacent_to_the_leap() {
    use siderust::formats::ilrs::cpf::{parse_cpf, write_cpf};
    for (mjd, sod, label) in [
        (57753, "86399.999999999", "2016-12-31T23:59:59.999999999Z"),
        (57754, "0.000000000", "2017-01-01T00:00:00Z"),
    ] {
        let input = format!("H1 CPF 2 TEST 2016 12 31 00 1\nH2 TEST 1 5 2016 12 31 00 00 00 2017 01 01 00 00 00 60 0 0\n10 1 {mjd} {sod} 0 1000000 2000000 3000000\n99\n");
        let cpf = parse_cpf(&input).unwrap();
        assert_eq!(
            cpf.ephemeris.entries[0].epoch.try_to_chrono().unwrap(),
            label.parse::<DateTime<Utc>>().unwrap()
        );
        let mut output = Vec::new();
        write_cpf(&cpf, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains(&format!("10 1 {mjd} {sod}")), "{output}");
        let roundtrip = parse_cpf(&output).unwrap();
        assert_eq!(
            cpf.ephemeris.entries[0].epoch,
            roundtrip.ephemeris.entries[0].epoch
        );
    }
}
