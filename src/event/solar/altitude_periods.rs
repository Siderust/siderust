// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! ## Scientific scope
//!
//! Sun‑specific routines for finding time intervals where the topocentric
//! altitude of the Sun is above, below, or within a given band — the
//! kinematic basis for day/night classification, twilight detection, and
//! observation‑window planning. The Sun position comes from
//! `Sun::get_horizontal` (VSOP87 + nutation + aberration); refraction is
//! not applied.
//!
//! ## Technical scope
//!
//! Crate‑internal API: `sun_altitude_rad`, [`find_day_periods`],
//! [`find_night_periods`], [`find_sun_range_periods`]. Period finding uses
//! Chebyshev-first labelled crossing discovery with precise validation and
//! local scan+Brent fallback via [`crate::event::search::crossings`].
//! Below‑threshold and range variants are derived via interval algebra.
//!
//! ## References
//! None.

use crate::bodies::solar_system::Sun;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::SearchOptsV2;
use crate::event::altitude::{CrossingDirection, CrossingEvent};
use crate::event::search::crossings;
use crate::event::search::intervals;
use crate::qtty::*;
use crate::time::{complement_within, Interval, JulianDate, ModifiedJulianDate};

// =============================================================================
// Constants
// =============================================================================

/// Scan step used only for scan+Brent fallback segments (2 hours in days).
const SCAN_STEP: Days = Quantity::<Hour>::new(2.0).to_const::<Day>();

// =============================================================================
// Core Altitude Function
// =============================================================================

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
/// Positive above the horizon, negative below.
pub(crate) fn sun_altitude_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Quantity<Radian> {
    let jd: JulianDate = mjd.to::<crate::JD>();
    Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .alt()
        .to::<Radian>()
}

// =============================================================================
// Main API
// =============================================================================

/// Finds day periods (Sun **above** `threshold`) inside `period`.
pub(crate) fn find_day_periods_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let thr = threshold.to::<Radian>();
    let signal = |t: ModifiedJulianDate| -> f64 { sun_altitude_rad(t, &site).sin() };
    let (labeled, start_above, _) =
        crossings::find_labelled_crossings(period, SCAN_STEP, &signal, thr.sin(), opts);
    intervals::build_above_periods_directed(&labeled, period, start_above)
}

/// Finds night periods (Sun **below** `twilight`) inside `period`.
pub(crate) fn find_night_periods_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    twilight: Degrees,
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let days = find_day_periods_with_search_opts(site, period, twilight, opts);
    complement_within(period, &days)
}

/// Finds periods where Sun altitude is within `range` `[min, max]`.
pub(crate) fn find_sun_range_periods_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    range: (Degrees, Degrees),
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let above_min = find_day_periods_with_search_opts(site, period, range.0, opts);
    let above_max = find_day_periods_with_search_opts(site, period, range.1, opts);
    let below_max = intervals::complement(period, &above_max);
    intervals::intersect(&above_min, &below_max)
}

pub(crate) fn find_sun_crossings_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
) -> Vec<CrossingEvent> {
    let thr = threshold.to::<Radian>();
    let signal = |t: ModifiedJulianDate| -> f64 { sun_altitude_rad(t, &site).sin() };
    let (labeled, _, _) =
        crossings::find_labelled_crossings(period, SCAN_STEP, &signal, thr.sin(), opts);
    labeled
        .iter()
        .map(|crossing| CrossingEvent {
            mjd: crossing.t,
            direction: if crossing.direction > 0 {
                CrossingDirection::Rising
            } else {
                CrossingDirection::Setting
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn generic_day_periods(
        site: Geodetic<ECEF>,
        period: Interval<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<Interval<ModifiedJulianDate>> {
        let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };
        intervals::above_threshold_periods(period, SCAN_STEP, &f, threshold.to::<Radian>())
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
    fn test_find_night_periods() {
        use crate::event::solar::twilight;

        let site = greenwich_site();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);
        let period = Interval::new(mjd_start, mjd_end);

        let nights = find_night_periods_with_search_opts(
            site,
            period,
            twilight::ASTRONOMICAL,
            SearchOptsV2::default(),
        );
        assert!(
            !nights.is_empty(),
            "Should find night periods at 51° latitude"
        );

        for night in &nights {
            assert!(
                (night.end.raw() - night.start.raw()) > Days::new(0.0),
                "Night duration should be positive"
            );
            assert!(
                (night.end.raw() - night.start.raw()) < Days::new(1.0),
                "Night should be less than 24 hours"
            );
        }
    }

    #[test]
    fn test_find_altitude_range_periods() {
        let site = greenwich_site();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);

        let period = Interval::new(mjd_start, mjd_end);

        let nights = find_sun_range_periods_with_search_opts(
            site,
            period,
            (Degrees::new(-90.0), Degrees::new(-18.0)),
            SearchOptsV2::default(),
        );

        assert!(!nights.is_empty(), "Should find night periods using range");

        let nautical = find_sun_range_periods_with_search_opts(
            site,
            period,
            (Degrees::new(-18.0), Degrees::new(-12.0)),
            SearchOptsV2::default(),
        );

        assert!(
            !nautical.is_empty(),
            "Should find nautical twilight periods"
        );
    }

    #[test]
    fn chebyshev_solar_threshold_periods_match_scan_baseline() {
        let site = greenwich_site();
        let period = Interval::new(utc_mjd(2026, 1, 1), utc_mjd(2026, 1, 8));

        for threshold in [
            Degrees::new(0.0),
            Degrees::new(-6.0),
            Degrees::new(-12.0),
            Degrees::new(-18.0),
        ] {
            let chebyshev_days =
                find_day_periods_with_search_opts(site, period, threshold, SearchOptsV2::default());
            let scan_days = generic_day_periods(site, period, threshold);
            assert_periods_close(&chebyshev_days, &scan_days);

            let chebyshev_nights = find_night_periods_with_search_opts(
                site,
                period,
                threshold,
                SearchOptsV2::default(),
            );
            let scan_nights = complement_within(period, &scan_days);
            assert_periods_close(&chebyshev_nights, &scan_nights);
        }
    }

    #[test]
    fn chebyshev_solar_range_periods_match_scan_baseline() {
        let site = greenwich_site();
        let period = Interval::new(utc_mjd(2026, 1, 1), utc_mjd(2026, 1, 8));
        let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };

        let chebyshev = find_sun_range_periods_with_search_opts(
            site,
            period,
            (Degrees::new(-18.0), Degrees::new(-12.0)),
            SearchOptsV2::default(),
        );
        let scan = intervals::in_range_periods(
            period,
            SCAN_STEP,
            &f,
            Degrees::new(-18.0).to::<Radian>(),
            Degrees::new(-12.0).to::<Radian>(),
        );

        assert_periods_close(&chebyshev, &scan);
    }

    #[test]
    fn chebyshev_solar_no_crossing_matches_scan_baseline() {
        let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(80.0), Meters::new(0.0));
        let period = Interval::new(utc_mjd(2026, 6, 20), utc_mjd(2026, 6, 27));
        let threshold = Degrees::new(-18.0);

        let chebyshev_nights =
            find_night_periods_with_search_opts(site, period, threshold, SearchOptsV2::default());
        let scan_days = generic_day_periods(site, period, threshold);
        let scan_nights = complement_within(period, &scan_days);

        assert_periods_close(&chebyshev_nights, &scan_nights);
        assert!(
            chebyshev_nights.is_empty(),
            "80°N near solstice should stay above -18°"
        );
    }
}
