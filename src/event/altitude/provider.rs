// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Provider, Trait-Based Dispatch Layer
//!
//! ## Scientific scope
//!
//! Provides a uniform, type-safe way to evaluate topocentric altitude *h(t)* for
//! the Sun, Moon, fixed stars, planets, and ICRS directions. Period semantics
//! (`above_threshold`, `below_threshold`, `altitude_ranges`) live in
//! [`super::events`]; this module routes single-point altitude evaluation and
//! optional internal search optimizations. Refraction is not applied.
//!
//! ## Technical scope
//!
//! Defines the [`AltitudeProvider`] trait and implementations that normalise
//! altitude evaluation across celestial targets. Optimized Sun/Moon period
//! search hooks are internal implementation details consumed by the event
//! engine in [`super::events`].

use super::search::{SearchOpts, SearchOptsV2};
use super::types::CrossingEvent;
use crate::bodies::solar_system;
use crate::bodies::Star;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::spherical::direction;
use crate::qtty::*;
use crate::time::{Interval, ModifiedJulianDate};

// Imports for planet altitude support
use crate::coordinates::{cartesian, centers::Geocentric, frames};
use crate::event::horizontal;
use crate::time::JulianDate;

// ---------------------------------------------------------------------------
// Trait Definition
// ---------------------------------------------------------------------------

/// Unified interface for evaluating topocentric altitude of any celestial target.
///
/// Implementors delegate single-point altitude to the appropriate analytical or
/// numerical engine. Period queries use the free functions
/// [`super::above_threshold`], [`super::below_threshold`], and
/// [`super::altitude_ranges`], which are generic over this trait.
///
/// Time scale note: all `ModifiedJulianDate` values are on the canonical JD(TT)
/// axis (`tempoch` semantics). Convert UTC instants with
/// `tempoch::Time::<tempoch::UTC>::from_chrono(...).to::<tempoch::TT>().into()`
/// into `ModifiedJulianDate` before using this API.
pub trait AltitudeProvider {
    /// Compute the altitude of this body at a single instant.
    ///
    /// Returns the topocentric altitude in radians.
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians;

    /// Compute the altitude with an explicit apparent-position correction policy.
    fn altitude_at_with_policy(
        &self,
        observer: &Geodetic<ECEF>,
        mjd: ModifiedJulianDate,
        policy: crate::astro::apparent::CorrectionPolicy,
    ) -> Radians {
        let _ = policy;
        self.altitude_at(observer, mjd)
    }

    /// Hint for the scan step to use when searching for events.
    ///
    /// Returns `None` to use the default (10 minutes). Bodies with slower
    /// apparent motion (like the Moon) can return a larger step for efficiency.
    fn scan_step_hint(&self) -> Option<Days> {
        None
    }

    /// Internal optimized above-threshold path with extended search options.
    fn above_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = (observer, window, threshold, opts);
        None
    }

    /// Internal optimized below-threshold path with extended search options.
    fn below_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = (observer, window, threshold, opts);
        None
    }

    /// Internal optimized altitude-range path with extended search options.
    fn altitude_range_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        min_altitude: Degrees,
        max_altitude: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = (observer, window, min_altitude, max_altitude, opts);
        None
    }

    /// Internal optimized threshold-crossing path with extended search options.
    fn crossings_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<CrossingEvent>> {
        let _ = (observer, window, threshold, opts);
        None
    }
}

// ---------------------------------------------------------------------------
// Implementations
// ---------------------------------------------------------------------------

/// **Sun** — specialized period search hooks delegate to [`crate::event::solar`].
impl AltitudeProvider for solar_system::Sun {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::event::solar::sun_altitude_rad(mjd, observer)
    }

    fn above_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(crate::event::solar::find_day_periods_with_search_opts(
            observer,
            window,
            threshold,
            SearchOptsV2::from_legacy(opts),
        ))
    }

    fn below_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(crate::event::solar::find_night_periods_with_search_opts(
            observer,
            window,
            threshold,
            SearchOptsV2::from_legacy(opts),
        ))
    }

    fn altitude_range_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        min_altitude: Degrees,
        max_altitude: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(
            crate::event::solar::find_sun_range_periods_with_search_opts(
                observer,
                window,
                (min_altitude, max_altitude),
                SearchOptsV2::from_legacy(opts),
            ),
        )
    }

    fn crossings_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<CrossingEvent>> {
        Some(crate::event::solar::find_sun_crossings_with_search_opts(
            observer,
            window,
            threshold,
            SearchOptsV2::from_legacy(opts),
        ))
    }
}

/// **Moon** — specialized period search hooks delegate to [`crate::event::lunar`].
impl AltitudeProvider for solar_system::Moon {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::event::lunar::moon_altitude_rad(mjd, observer)
    }

    fn above_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(
            crate::event::lunar::find_moon_above_horizon_with_search_opts(
                observer,
                window,
                threshold,
                SearchOptsV2::from_legacy(opts),
            ),
        )
    }

    fn below_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(
            crate::event::lunar::find_moon_below_horizon_with_search_opts(
                observer,
                window,
                threshold,
                SearchOptsV2::from_legacy(opts),
            ),
        )
    }

    fn altitude_range_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        min_altitude: Degrees,
        max_altitude: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        Some(
            crate::event::lunar::find_moon_altitude_range_with_search_opts(
                observer,
                window,
                (min_altitude, max_altitude),
                SearchOptsV2::from_legacy(opts),
            ),
        )
    }

    fn crossings_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<CrossingEvent>> {
        Some(crate::event::lunar::find_moon_crossings_with_search_opts(
            observer,
            window,
            threshold,
            SearchOptsV2::from_legacy(opts),
        ))
    }

    fn scan_step_hint(&self) -> Option<Days> {
        // Moon moves slower, 2-hour steps are sufficient
        Some(Hours::new(2.0).to::<Day>())
    }
}

/// **Star** — delegates single-point altitude to [`direction::ICRS`].
impl AltitudeProvider for Star<'_> {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        let dir = direction::ICRS::from(self);
        dir.altitude_at(observer, mjd)
    }

    fn altitude_at_with_policy(
        &self,
        observer: &Geodetic<ECEF>,
        mjd: ModifiedJulianDate,
        policy: crate::astro::apparent::CorrectionPolicy,
    ) -> Radians {
        let dir = direction::ICRS::from(self);
        dir.altitude_at_with_policy(observer, mjd, policy)
    }

    fn above_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = opts;
        let dir = direction::ICRS::from(self);
        dir.above_threshold_search(observer, window, threshold, opts)
    }

    fn below_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let dir = direction::ICRS::from(self);
        dir.below_threshold_search(observer, window, threshold, opts)
    }

    fn altitude_range_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        min_altitude: Degrees,
        max_altitude: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let dir = direction::ICRS::from(self);
        dir.altitude_range_search(observer, window, min_altitude, max_altitude, opts)
    }
}

/// **direction::ICRS** — fixed-sky direction altitude evaluation.
impl AltitudeProvider for direction::ICRS {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::event::stellar::fixed_star_altitude_rad(mjd, observer, self.ra(), self.dec())
    }

    fn altitude_at_with_policy(
        &self,
        observer: &Geodetic<ECEF>,
        mjd: ModifiedJulianDate,
        policy: crate::astro::apparent::CorrectionPolicy,
    ) -> Radians {
        crate::event::stellar::fixed_star_altitude_rad_with_policy(
            mjd,
            observer,
            self.ra(),
            self.dec(),
            policy,
        )
    }

    fn above_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = opts;
        Some(crate::event::stellar::find_star_above_periods(
            self.ra(),
            self.dec(),
            observer,
            window,
            threshold,
        ))
    }

    fn below_threshold_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        threshold: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = opts;
        Some(crate::event::stellar::find_star_below_periods(
            self.ra(),
            self.dec(),
            observer,
            window,
            threshold,
        ))
    }

    fn altitude_range_search(
        &self,
        observer: Geodetic<ECEF>,
        window: Interval<ModifiedJulianDate>,
        min_altitude: Degrees,
        max_altitude: Degrees,
        opts: SearchOpts,
    ) -> Option<Vec<Interval<ModifiedJulianDate>>> {
        let _ = opts;
        Some(crate::event::stellar::find_star_range_periods(
            self.ra(),
            self.dec(),
            observer,
            window,
            (min_altitude, max_altitude),
        ))
    }
}

// ---------------------------------------------------------------------------
// Implementations: VSOP87 Planets (Mercury–Neptune)
// ---------------------------------------------------------------------------

use crate::coordinates::transform::Transform;

/// Scan step for planet altitude threshold detection (2 hours in days).
///
/// Planets' apparent motion is dominated by Earth's rotation, so the same
/// 2‑hour scan step used for the Sun is adequate.
const PLANET_SCAN_STEP: Days = Quantity::<Hour>::new(2.0).to_const::<Day>();

/// Computes the topocentric altitude (in radians) of a VSOP87 planet at a
/// given instant, using the full VSOP87 → geocentric equatorial → topocentric
/// → horizontal pipeline.
fn vsop87_planet_altitude_rad<F>(
    vsop87e_fn: F,
    mjd: ModifiedJulianDate,
    site: &Geodetic<ECEF>,
) -> Radians
where
    F: Fn(
        JulianDate,
    ) -> cartesian::Position<
        crate::coordinates::centers::Barycentric,
        frames::EclipticMeanJ2000,
        AstronomicalUnit,
    >,
{
    let jd: JulianDate = mjd.to::<crate::JD>();
    // 1) VSOP87e → barycentric ecliptic J2000
    let bary_ecl = vsop87e_fn(jd);
    // 2) Frame rotation + center shift → geocentric equatorial J2000
    let geo_equ: cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, AstronomicalUnit> =
        bary_ecl.transform(jd);
    // 3–4) Topocentric parallax + precession/nutation → true‑of‑date RA/Dec,
    //       then equatorial → horizontal
    let topo = horizontal::geocentric_j2000_to_apparent_topocentric(&geo_equ, *site, jd);
    let horiz = horizontal::equatorial_to_horizontal(&topo, *site, jd);
    horiz.alt().to::<Radian>()
}

/// Helper macro: implement [`AltitudeProvider`] for a VSOP87‑backed planet.
macro_rules! impl_altitude_provider_vsop87 {
    ($($Planet:ident),+ $(,)?) => {
        $(
            impl AltitudeProvider for solar_system::$Planet {
                fn altitude_at(
                    &self,
                    observer: &Geodetic<ECEF>,
                    mjd: ModifiedJulianDate,
                ) -> Radians {
                    vsop87_planet_altitude_rad(
                        solar_system::$Planet::vsop87e, mjd, observer,
                    )
                }

                fn scan_step_hint(&self) -> Option<Days> {
                    Some(PLANET_SCAN_STEP)
                }
            }
        )+
    };
}

impl_altitude_provider_vsop87!(Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog;
    use crate::event::altitude::{above_threshold, altitude_ranges, below_threshold, SearchOpts};

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn one_day_window() -> Interval<ModifiedJulianDate> {
        Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60001.0),
        )
    }

    fn one_week_window() -> Interval<ModifiedJulianDate> {
        Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60007.0),
        )
    }

    #[test]
    fn sun_above_horizon_via_standard_api() {
        let periods = above_threshold(
            &solar_system::Sun,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(!periods.is_empty(), "Sun should be above horizon at 51°N");
        for p in &periods {
            assert!(p.length() > Days::new(0.0));
            assert!(p.length() < Days::new(1.0));
        }
    }

    #[test]
    fn moon_above_horizon_via_standard_api() {
        let periods = above_threshold(
            &solar_system::Moon,
            &greenwich(),
            one_week_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(
            !periods.is_empty(),
            "Moon should be above horizon at some point in a week"
        );
    }

    #[test]
    fn star_above_horizon_via_standard_api() {
        let sirius = &catalog::SIRIUS;
        let periods = above_threshold(
            sirius,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(
            !periods.is_empty(),
            "Sirius should be above horizon for part of the day"
        );
    }

    #[test]
    fn icrs_direction_above_horizon_via_standard_api() {
        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let periods = above_threshold(
            &sirius_dir,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(
            !periods.is_empty(),
            "direction::ICRS for Sirius should match Star result"
        );
    }

    #[test]
    fn star_and_icrs_direction_agree() {
        let sirius = &catalog::SIRIUS;
        let sirius_dir = direction::ICRS::from(sirius);
        let window = one_day_window();
        let observer = greenwich();
        let opts = SearchOpts::default();

        let star_periods = above_threshold(sirius, &observer, window, Degrees::new(0.0), opts);
        let dir_periods = above_threshold(&sirius_dir, &observer, window, Degrees::new(0.0), opts);

        assert_eq!(
            star_periods.len(),
            dir_periods.len(),
            "Star and direction::ICRS should produce the same number of periods"
        );
        for (sp, dp) in star_periods.iter().zip(dir_periods.iter()) {
            assert!(
                (sp.start.raw() - dp.start.raw()).abs() < Days::new(1e-6),
                "Period starts should match"
            );
            assert!(
                (sp.end.raw() - dp.end.raw()).abs() < Days::new(1e-6),
                "Period ends should match"
            );
        }
    }

    #[test]
    fn altitude_at_consistent_across_types() {
        let observer = greenwich();
        let mjd = crate::time::ModifiedJulianDate::new(51544.5);

        let sun_alt = solar_system::Sun.altitude_at(&observer, mjd);
        assert!(sun_alt.abs() < Radians::new(std::f64::consts::FRAC_PI_2));

        let moon_alt = solar_system::Moon.altitude_at(&observer, mjd);
        assert!(moon_alt.abs() < Radians::new(std::f64::consts::FRAC_PI_2));

        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let star_alt = sirius_dir.altitude_at(&observer, mjd);
        assert!(star_alt.abs() < Radians::new(std::f64::consts::FRAC_PI_2));
    }

    #[test]
    fn full_sky_range_returns_full_window() {
        let periods = altitude_ranges(
            &solar_system::Sun,
            &greenwich(),
            one_day_window(),
            Degrees::new(-90.0),
            Degrees::new(90.0),
            SearchOpts::default(),
        );
        assert!(
            !periods.is_empty(),
            "Full sky range should return at least one period"
        );
        let total: f64 = periods.iter().map(|p| p.length().value()).sum();
        assert!(
            (total - 1.0).abs() < 0.01,
            "Full sky range should span ~1 day, got {} days",
            total
        );
    }

    #[test]
    fn polaris_circumpolar_via_standard_api() {
        let polaris = &catalog::POLARIS;
        let periods = above_threshold(
            polaris,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert_eq!(
            periods.len(),
            1,
            "Polaris should be continuously above horizon at 51°N"
        );
        assert!(
            ((periods[0].end.raw() - periods[0].start.raw()) - Days::new(1.0)).abs()
                < Days::new(0.01),
            "Polaris up-period should span the full day"
        );
    }

    #[test]
    fn polaris_never_below_minus80_via_standard_api() {
        let polaris = &catalog::POLARIS;
        let periods = below_threshold(
            polaris,
            &greenwich(),
            one_day_window(),
            Degrees::new(-80.0),
            SearchOpts::default(),
        );
        assert!(
            periods.is_empty(),
            "Polaris should never be below -80° at 51°N"
        );
    }

    #[test]
    fn empty_window_returns_empty() {
        let window = Interval::new(
            crate::time::ModifiedJulianDate::new(60000.0),
            crate::time::ModifiedJulianDate::new(60000.0),
        );
        let periods = above_threshold(
            &solar_system::Sun,
            &greenwich(),
            window,
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(periods.is_empty(), "Empty window should return no periods");
    }

    #[test]
    fn below_threshold_sun_night_via_standard_api() {
        let nights = below_threshold(
            &solar_system::Sun,
            &greenwich(),
            one_week_window(),
            Degrees::new(-18.0),
            SearchOpts::default(),
        );
        assert!(!nights.is_empty(), "Should find astronomical night at 51°N");
    }

    #[test]
    fn altitude_range_twilight_via_standard_api() {
        let bands = altitude_ranges(
            &solar_system::Sun,
            &greenwich(),
            Interval::new(
                crate::time::ModifiedJulianDate::new(60000.0),
                crate::time::ModifiedJulianDate::new(60002.0),
            ),
            Degrees::new(-18.0),
            Degrees::new(-12.0),
            SearchOpts::default(),
        );
        assert!(
            bands.len() >= 2,
            "Should find at least 2 twilight bands in 2 days, found {}",
            bands.len()
        );
    }

    #[test]
    fn periods_are_sorted_and_non_overlapping() {
        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let periods = above_threshold(
            &sirius_dir,
            &greenwich(),
            one_week_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        for w in periods.windows(2) {
            assert!(
                w[0].end <= w[1].start,
                "Periods should be non-overlapping and sorted: {:?} vs {:?}",
                w[0],
                w[1]
            );
        }
    }

    // --- Planet altitude ---

    #[test]
    fn mars_altitude_at_is_finite() {
        let alt = solar_system::Mars
            .altitude_at(&greenwich(), crate::time::ModifiedJulianDate::new(60000.5));
        assert!(alt.is_finite());
        assert!(
            alt.abs() < Radians::new(std::f64::consts::FRAC_PI_2),
            "Mars altitude should be within ±90°"
        );
    }

    #[test]
    fn jupiter_above_horizon_via_standard_api() {
        let periods = above_threshold(
            &solar_system::Jupiter,
            &greenwich(),
            one_week_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(
            !periods.is_empty(),
            "Jupiter should be above horizon at some point in a week at 51°N"
        );
    }

    #[test]
    fn planet_altitudes_are_realistic() {
        let observer = greenwich();
        let mjd = crate::time::ModifiedJulianDate::new(60000.5);
        // All planets should return finite altitudes
        let mercury_alt = solar_system::Mercury.altitude_at(&observer, mjd);
        let venus_alt = solar_system::Venus.altitude_at(&observer, mjd);
        let mars_alt = solar_system::Mars.altitude_at(&observer, mjd);
        let saturn_alt = solar_system::Saturn.altitude_at(&observer, mjd);
        let uranus_alt = solar_system::Uranus.altitude_at(&observer, mjd);
        let neptune_alt = solar_system::Neptune.altitude_at(&observer, mjd);

        for (name, alt) in [
            ("Mercury", mercury_alt),
            ("Venus", venus_alt),
            ("Mars", mars_alt),
            ("Saturn", saturn_alt),
            ("Uranus", uranus_alt),
            ("Neptune", neptune_alt),
        ] {
            assert!(alt.is_finite(), "{name} altitude should be finite");
            assert!(
                alt.abs() < Radians::new(std::f64::consts::FRAC_PI_2),
                "{name} altitude should be within ±90°, got {alt}"
            );
        }
    }
}
