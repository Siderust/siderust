// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Scales, Formats, and Period Conversions Example
//!
//! Run with: `cargo run --example 05_time_periods`
//!
//! This example demonstrates `tempoch` (re-exported as `siderust::time`):
//! - Constructing an instant from `chrono::DateTime<Utc>`
//! - Viewing the same absolute instant in each supported time scale:
//!   `JD`, `JDE`, `MJD`, `TDB`, `TT`, `TAI`, `TCG`, `TCB`, `GPS`, `UnixTime`, `UT`
//! - Using the common type aliases: `JulianDate`, `JulianEphemerisDay`,
//!   `ModifiedJulianDate`, `UniversalTime`
//! - Converting `Period<S>` between scales and to/from `UtcPeriod`
//!
//! Notes:
//! - Scale conversions route through the canonical `JD(TT)` representation.
//! - `UT` is Earth-rotation based; `Time::<UT>::delta_t()` exposes `ΔT = TT − UT`.

use chrono::{DateTime, Duration, Utc};
use qtty::{Days, Second};
use siderust::time::{
    Interval, JulianDate, JulianEphemerisDay, ModifiedJulianDate, Period, Time, TimeScale,
    UniversalTime, UnixTime, UtcPeriod, GPS, JD, JDE, MJD, TAI, TCB, TCG, TDB, TT, UT,
};

fn print_scale<S: TimeScale>(label: &str, time: Time<S>, reference_jd: JulianDate) {
    let jd_back: JulianDate = time.to::<JD>();
    let drift_s = (jd_back - reference_jd).to::<Second>();
    println!(
        "   {:<8} value = {:>16.9}  | JD roundtrip drift = {:>11.3e} s",
        label, time, drift_s
    );
}

fn print_period<S: TimeScale>(label: &str, period: &Period<S>) {
    println!(
        "   {:<8} [{:>16.9}]  Δ = {}",
        label,
        period,
        period.duration_days()
    );
}

fn main() {
    println!("Time Scales, Formats, and Period Conversions");
    println!("============================================\n");

    let utc_ref = DateTime::<Utc>::from_timestamp(946_728_000, 0).expect("valid UTC timestamp");
    let jd: JulianDate = JulianDate::from_utc(utc_ref);

    let jde: JulianEphemerisDay = jd.to::<JDE>();
    let mjd: ModifiedJulianDate = jd.to::<MJD>();
    let tdb: Time<TDB> = jd.to::<TDB>();
    let tt: Time<TT> = jd.to::<TT>();
    let tai: Time<TAI> = jd.to::<TAI>();
    let tcg: Time<TCG> = jd.to::<TCG>();
    let tcb: Time<TCB> = jd.to::<TCB>();
    let gps: Time<GPS> = jd.to::<GPS>();
    let unix: Time<UnixTime> = jd.to::<UnixTime>();
    let ut: UniversalTime = jd.to::<UT>();

    println!("Reference UTC instant: {}\n", utc_ref.to_rfc3339());

    println!("1) Each supported time scale for the same instant:");
    print_scale("JD", jd, jd);
    print_scale("JDE", jde, jd);
    print_scale("MJD", mjd, jd);
    print_scale("TDB", tdb, jd);
    print_scale("TT", tt, jd);
    print_scale("TAI", tai, jd);
    print_scale("TCG", tcg, jd);
    print_scale("TCB", tcb, jd);
    print_scale("GPS", gps, jd);
    print_scale("Unix", unix, jd);
    print_scale("UT", ut, jd);
    println!(
        "   {:<8} delta_t = {:.3} s (TT - UT)\n",
        "UT",
        ut.delta_t().value()
    );

    println!("2) Time formats / aliases:");
    println!("   JulianDate alias:         {}", jd);
    println!("   JulianEphemerisDay alias: {}", jde);
    println!("   ModifiedJulianDate alias: {}", mjd);
    println!("   UniversalTime alias:      {}", ut);
    let utc_roundtrip = jd.to_utc().expect("JD should convert back to UTC");
    println!(
        "   UTC roundtrip from JD:    {}\n",
        utc_roundtrip.to_rfc3339()
    );

    println!("3) Period representations and conversions:");
    let period_jd: Period<JD> = Period::new(jd, jd + Days::new(0.5));
    let period_jde: Period<JDE> = period_jd.to::<JDE>().expect("JD -> JDE period conversion");
    let period_mjd: Period<MJD> = period_jd.to::<MJD>().expect("JD -> MJD period conversion");
    let period_tdb: Period<TDB> = period_jd.to::<TDB>().expect("JD -> TDB period conversion");
    let period_tt: Period<TT> = period_jd.to::<TT>().expect("JD -> TT period conversion");
    let period_tai: Period<TAI> = period_jd.to::<TAI>().expect("JD -> TAI period conversion");
    let period_tcg: Period<TCG> = period_jd.to::<TCG>().expect("JD -> TCG period conversion");
    let period_tcb: Period<TCB> = period_jd.to::<TCB>().expect("JD -> TCB period conversion");
    let period_gps: Period<GPS> = period_jd.to::<GPS>().expect("JD -> GPS period conversion");
    let period_unix: Period<UnixTime> = period_jd
        .to::<UnixTime>()
        .expect("JD -> Unix period conversion");
    let period_ut: Period<UT> = period_jd.to::<UT>().expect("JD -> UT period conversion");
    let period_utc: UtcPeriod = period_jd
        .to::<DateTime<Utc>>()
        .expect("JD -> UTC period conversion");

    print_period("JD", &period_jd);
    print_period("JDE", &period_jde);
    print_period("MJD", &period_mjd);
    print_period("TDB", &period_tdb);
    print_period("TT", &period_tt);
    print_period("TAI", &period_tai);
    print_period("TCG", &period_tcg);
    print_period("TCB", &period_tcb);
    print_period("GPS", &period_gps);
    print_period("Unix", &period_unix);
    print_period("UT", &period_ut);
    println!(
        "   {:<8} [{} -> {}]  Δ = {:.6} days ({:.0} s)\n",
        "UTC",
        period_utc.start.to_rfc3339(),
        period_utc.end.to_rfc3339(),
        period_utc.duration_days(),
        period_utc.duration_seconds()
    );

    println!("4) UtcPeriod / Interval<DateTime<Utc>> conversions back to typed periods:");
    let utc_window: UtcPeriod = Interval::new(utc_ref, utc_ref + Duration::hours(6));
    let from_utc_jd: Period<JD> = utc_window.to::<JD>();
    let from_utc_mjd: Period<MJD> = utc_window.to::<MJD>();
    let from_utc_ut: Period<UT> = utc_window.to::<UT>();
    let from_utc_unix: Period<UnixTime> = utc_window.to::<UnixTime>();

    println!(
        "   UTC      [{} -> {}]  Δ = {:.6} days",
        utc_window.start.to_rfc3339(),
        utc_window.end.to_rfc3339(),
        utc_window.duration_days()
    );
    print_period("JD", &from_utc_jd);
    print_period("MJD", &from_utc_mjd);
    print_period("UT", &from_utc_ut);
    print_period("Unix", &from_utc_unix);

    let utc_roundtrip_period: UtcPeriod = from_utc_mjd
        .to::<DateTime<Utc>>()
        .expect("MJD -> UTC period conversion");
    println!(
        "   UTC<-MJD [{} -> {}]",
        utc_roundtrip_period.start.to_rfc3339(),
        utc_roundtrip_period.end.to_rfc3339()
    );
}
