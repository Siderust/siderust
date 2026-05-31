// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(missing_docs)]

use qtty::time::JULIAN_YEAR;
use siderust::qtty::Days;
use siderust::time::JulianDate;

#[test]
fn julian_date_arithmetic_display_and_chrono_roundtrip() {
    let mut jd = JulianDate::try_new(Days::new(2_450_000.0)).unwrap();
    let printed = format!("{jd}");
    assert_eq!(printed, format!("{}", jd.raw()));

    jd += Days::new(2.0);
    jd -= Days::new(0.5);

    let julian_year_days = JULIAN_YEAR.to::<qtty::unit::Day>();
    let with_years = jd + julian_year_days;
    let day_span: Days = with_years.raw() - jd.raw();
    assert!((day_span.value() - julian_year_days.value()).abs() < 1e-9);

    let min = if with_years < jd { with_years } else { jd };
    assert_eq!(min, jd);

    let utc = jd
        .to::<siderust::time::UTC>()
        .to_chrono()
        .expect("valid UTC");
    let roundtrip: JulianDate = siderust::time::JulianDate::from(utc);
    assert!((roundtrip.raw().value() - jd.raw().value()).abs() < 1e-6);
}
