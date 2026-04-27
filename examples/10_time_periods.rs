// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time wrappers and interval helpers.

use chrono::{Duration, Utc};
use siderust::qtty::Days;
use siderust::time::{Interval, JulianDate, ModifiedJulianDate, UTC};

fn main() {
    let now_utc = Utc::now();
    let jd = JulianDate::from_utc(now_utc);
    let mjd: ModifiedJulianDate = jd.into();

    println!("UTC now : {}", now_utc.to_rfc3339());
    println!("JD (TT) : {}", jd);
    println!("MJD(TT) : {}", mjd);
    println!("Back UTC: {}", jd.to_utc().to_rfc3339());

    let tomorrow = jd + Days::new(1.0);
    let window = Interval::<JulianDate>::new(jd, tomorrow);
    println!("1-day window length: {}", window.end - window.start);

    let utc_time = tempoch::Time::<UTC>::from_chrono(now_utc);
    let utc_window = Interval::<tempoch::Time<UTC>>::new(utc_time, tempoch::Time::<UTC>::from_chrono(now_utc + Duration::hours(6)));
    println!("UTC scale window: {}", utc_window.end - utc_window.start);
}
