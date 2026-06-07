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
use crate::event::search::intervals;
use crate::qtty::*;
use crate::time::{complement_within, Interval, JulianDate, ModifiedJulianDate};

use super::daily_events::solar_daily_crossing_events_impl;
use super::daily_events::{solar_daily_crossings_for_thresholds_impl, solar_daily_crossings_impl};

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
pub(crate) fn sun_altitude_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Quantity<Radian> {
    let jd: JulianDate = mjd.to::<crate::JD>();
    crate::bodies::solar_system::Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .alt()
        .to::<Radian>()
}

pub(crate) fn solar_above_threshold_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let (labelled, start_above) = solar_labelled_crossings(site, period, threshold, opts);
    intervals::build_above_periods_directed(&labelled, period, start_above)
}

pub(crate) fn solar_below_threshold_impl(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    twilight: Degrees,
    opts: InternalSearchConfig,
) -> Vec<Interval<ModifiedJulianDate>> {
    let days = solar_above_threshold_impl(site, period, twilight, opts);
    complement_within(period, &days)
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

    let start_above_min = sun_altitude_rad(period.start, &site).sin() > min_sin;
    let start_above_max = sun_altitude_rad(period.start, &site).sin() > max_sin;

    let above_min = intervals::build_above_periods_directed(min_crossings, period, start_above_min);
    let above_max = intervals::build_above_periods_directed(max_crossings, period, start_above_max);
    let below_max = intervals::complement(period, &above_max);
    intervals::intersect(&above_min, &below_max)
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
) -> (Vec<intervals::LabeledCrossing>, bool) {
    let thr_sin = threshold.to::<Radian>().sin();
    let start_above = sun_altitude_rad(period.start, &site).sin() > thr_sin;
    let (labelled, _) = solar_daily_crossings_impl(site, period, threshold, opts);
    (labelled, start_above)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::event::altitude::search::DIURNAL_CHEBY_SCAN_STEP;
    use chrono::{TimeZone, Utc};

    fn greenwich_site() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    fn utc_mjd(year: i32, month: u32, day: u32) -> ModifiedJulianDate {
        ModifiedJulianDate::from(
            Utc.with_ymd_and_hms(year, month, day, 0, 0, 0)
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
        let period = Interval::new(utc_mjd(2026, 1, 1), utc_mjd(2026, 1, 8));
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
        let period = Interval::new(utc_mjd(2026, 1, 1), utc_mjd(2026, 1, 8));

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
}
