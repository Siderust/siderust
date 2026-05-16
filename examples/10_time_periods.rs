// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time wrappers and interval helpers.

use chrono::{Duration, Utc};
use siderust::qtty::Days;
use siderust::time::{Interval, JulianDate, ModifiedJulianDate, UTC};

fn main() {
    let now_utc = Utc::now();
    let jd: JulianDate = siderust::time::Time::<siderust::time::UTC>::from_chrono(now_utc)
        .to::<siderust::time::TT>()
        .into();
    let mjd: ModifiedJulianDate = jd.to::<siderust::time::MJD>();

    println!("UTC now : {}", now_utc.to_rfc3339());
    println!("JD (TT) : {}", jd);
    println!("MJD(TT) : {}", mjd);
    println!(
        "Back UTC: {}",
        jd.to::<siderust::time::UTC>()
            .to_chrono()
            .map(|dt| dt.to_rfc3339())
            .unwrap_or_else(|| "N/A".into())
    );

    let tomorrow = jd + Days::new(1.0);
    let window = Interval::<JulianDate>::new(jd, tomorrow);
    println!(
        "1-day window length: {}",
        window.end.raw() - window.start.raw()
    );

    let utc_time = tempoch::Time::<UTC>::from_chrono(now_utc);
    let utc_window = Interval::<tempoch::Time<UTC>>::new(
        utc_time,
        tempoch::Time::<UTC>::from_chrono(now_utc + Duration::hours(6)),
    );
    println!("UTC scale window: {}", utc_window.end - utc_window.start);
}
