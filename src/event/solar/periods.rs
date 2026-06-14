// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! Crate-internal API built on the solar daily predictor with local Chebyshev
//! and scan+Brent fallback.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::InternalSearchConfig;
use crate::event::altitude::CrossingEvent;
use crate::event::search::periods as threshold_periods;
use crate::qtty::*;
use crate::time::{Interval, ModifiedJulianDate};

use super::altitude::sun_altitude_rad;
use super::daily_events::solar_daily_crossing_events_impl;
use super::daily_events::{solar_daily_crossings_for_thresholds_impl, solar_daily_crossings_impl};

pub(crate) fn solar_above_threshold_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let (labelled, start_above) = solar_labelled_crossings(site, period, threshold, opts);
    threshold_periods::assemble_above_threshold_periods(&labelled, period, start_above)
}

pub(crate) fn solar_below_threshold_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    twilight: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let above = solar_above_threshold_impl(site, period, twilight, opts);
    threshold_periods::complement_threshold_periods(period, &above)
}

pub(crate) fn solar_altitude_ranges_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    range: (Degrees, Degrees),
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    if period.end <= period.start {
        return Vec::new();
    }

    let crossings =
        solar_daily_crossings_for_thresholds_impl(site, period, &[range.0, range.1], opts);
    let min_crossings = &crossings[0];
    let max_crossings = &crossings[1];

    let min_sin = range.0.to::<Radian>().sin();
    let max_sin = range.1.to::<Radian>().sin();

    let start_sin_altitude = sun_altitude_rad(period.start, &site).sin();
    threshold_periods::assemble_in_range_periods(
        min_crossings,
        start_sin_altitude > min_sin,
        max_crossings,
        start_sin_altitude > max_sin,
        period,
    )
}

pub(crate) fn solar_crossings_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<CrossingEvent> {
    solar_daily_crossing_events_impl(site, period, threshold, opts)
}

fn solar_labelled_crossings(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> (Vec<crate::event::search::intervals::LabeledCrossing>, bool) {
    let thr_sin = threshold.to::<Radian>().sin();
    let start_above = sun_altitude_rad(period.start, &site).sin() > thr_sin;
    let (labelled, _) = solar_daily_crossings_impl(site, period, threshold, opts);
    (labelled, start_above)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::event::altitude::search::DIURNAL_CHEBY_SCAN_STEP;
    use crate::event::search::intervals;
    use crate::time::complement_within;
    use chrono::{TimeZone, Utc};

    fn greenwich_site() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    fn cta_s_site() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(-70.406944),
            Degrees::new(-24.627222),
            Meters::new(2100.0),
        )
    }

    fn utc_midnight_as_tt_mjd(year: i32, month: u32, day: u32) -> ModifiedJulianDate {
        ModifiedJulianDate::from(
            Utc.with_ymd_and_hms(year, month, day, 0, 0, 0)
                .single()
                .unwrap(),
        )
    }

    fn utc_datetime_as_tt_mjd(
        year: i32,
        month: u32,
        day: u32,
        hour: u32,
        minute: u32,
        second: u32,
    ) -> ModifiedJulianDate {
        ModifiedJulianDate::from(
            Utc.with_ymd_and_hms(year, month, day, hour, minute, second)
                .single()
                .unwrap(),
        )
    }

    fn assert_periods_close(
        actual: &[Interval<ModifiedJulianDate>],
        expected: &[Interval<ModifiedJulianDate>],
    ) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "actual={actual:?} expected={expected:?}"
        );
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert!(
                (a.start.raw() - e.start.raw()).abs() < Days::new(1e-6),
                "start mismatch: actual={a:?} expected={e:?}"
            );
            assert!(
                (a.end.raw() - e.end.raw()).abs() < Days::new(1e-6),
                "end mismatch: actual={a:?} expected={e:?}"
            );
        }
    }

    fn scan_day_periods(
        site: Geodetic<ECEF>,
        period: Interval<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<Interval<ModifiedJulianDate>> {
        let scan_step = DIURNAL_CHEBY_SCAN_STEP;
        let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };
        intervals::above_threshold_periods(period, scan_step, &f, threshold.to::<Radian>())
    }

    fn independent_solar_altitude_range_reference(
        site: Geodetic<ECEF>,
        period: Interval<ModifiedJulianDate>,
        range: (Degrees, Degrees),
        opts: InternalSearchConfig,
    ) -> Vec<Interval<ModifiedJulianDate>> {
        let above_min = solar_above_threshold_impl(site, period, range.0, opts);
        let above_max = solar_above_threshold_impl(site, period, range.1, opts);
        let below_max = intervals::complement(period, &above_max);
        intervals::intersect(&above_min, &below_max)
    }

    #[test]
    fn solar_altitude_range_batch_matches_independent_thresholds() {
        let site = greenwich_site();
        let period = Interval::new(
            utc_midnight_as_tt_mjd(2026, 1, 1),
            utc_midnight_as_tt_mjd(2026, 1, 8),
        );
        let range = (Degrees::new(-18.0), Degrees::new(-12.0));
        let opts = InternalSearchConfig::default();

        let batch = solar_altitude_ranges_impl(site, period, range, opts);
        let reference = independent_solar_altitude_range_reference(site, period, range, opts);
        assert_periods_close(&batch, &reference);
    }

    #[test]
    fn test_sun_altitude_basic() {
        let site = greenwich_site();
        let mjd: ModifiedJulianDate = crate::J2000.to::<crate::MJD>();
        let alt = sun_altitude_rad(mjd, &site);
        assert!(
            alt > Radians::new(-std::f64::consts::FRAC_PI_2)
                && alt < Radians::new(std::f64::consts::FRAC_PI_2)
        );
    }

    #[test]
    fn test_solar_below_threshold_periods() {
        use crate::event::solar::twilight;

        let site = greenwich_site();
        let period = Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60007.0),
        );

        let nights = solar_below_threshold_impl(
            site,
            period,
            twilight::ASTRONOMICAL,
            InternalSearchConfig::default(),
        );
        assert!(!nights.is_empty());
    }

    #[test]
    fn daily_predictor_periods_match_scan_baseline() {
        let site = greenwich_site();
        let period = Interval::new(
            utc_midnight_as_tt_mjd(2026, 1, 1),
            utc_midnight_as_tt_mjd(2026, 1, 8),
        );

        for threshold in [
            Degrees::new(0.0),
            Degrees::new(-6.0),
            Degrees::new(-12.0),
            Degrees::new(-18.0),
        ] {
            let daily = solar_above_threshold_impl(
                site,
                period,
                threshold,
                InternalSearchConfig::default(),
            );
            let scan = scan_day_periods(site, period, threshold);
            assert_periods_close(&daily, &scan);

            let daily_nights = solar_below_threshold_impl(
                site,
                period,
                threshold,
                InternalSearchConfig::default(),
            );
            let scan_nights = complement_within(period, &scan);
            assert_periods_close(&daily_nights, &scan_nights);
        }
    }

    #[test]
    fn cta_s_astronomical_nights_match_scan_baseline_2025() {
        let site = cta_s_site();
        let start = utc_datetime_as_tt_mjd(2025, 1, 1, 12, 0, 0);
        let end = utc_datetime_as_tt_mjd(2025, 7, 31, 12, 0, 0);
        let window = Interval::new(start, end);
        let threshold = Degrees::new(-18.0);

        let daily_days =
            solar_above_threshold_impl(site, window, threshold, InternalSearchConfig::default());
        let scan_days = scan_day_periods(site, window, threshold);

        assert_periods_close(&daily_days, &scan_days);

        let daily_nights =
            solar_below_threshold_impl(site, window, threshold, InternalSearchConfig::default());
        let scan_nights = complement_within(window, &scan_days);

        assert_periods_close(&daily_nights, &scan_nights);
        assert!(
            daily_nights.len() > 180,
            "CTA-S should have many astronomical-night periods in Jan-Jul 2025, got {}",
            daily_nights.len()
        );

        for night in &daily_nights {
            assert!(
                night.end > night.start,
                "night interval must have positive duration: {night:?}"
            );
            assert!(
                night.start >= window.start,
                "night starts before query window: {night:?}"
            );
            assert!(
                night.end <= window.end,
                "night ends after query window: {night:?}"
            );
            let midpoint =
                ModifiedJulianDate::new((night.start.raw() + night.end.raw()).value() / 2.0);
            assert!(
                sun_altitude_rad(midpoint, &site) < threshold.to::<Radian>(),
                "night midpoint must be below astronomical threshold: {night:?}"
            );
        }
    }
}
