// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Periods — Trait‑Based Dispatch Layer
//!
//! This module defines the [`AltitudePeriodsProvider`] trait and implementations
//! that normalise the "altitude periods" API across Sun, Moon, and fixed stars.
//!
//! ## Design
//!
//! The trait layer is the *dispatch* layer — no astronomical math lives
//! here.  Each `impl` delegates to the appropriate engine inside
//! [`calculus::solar`], [`calculus::lunar`], or [`calculus::stellar`].
//!
//! | Body type | Engine |
//! |-----------|--------|
//! | [`solar_system::Sun`] | [`calculus::solar::find_day_periods`] / [`calculus::solar::find_sun_range_periods`] |
//! | [`solar_system::Moon`] | [`calculus::lunar::find_moon_above_horizon`] / [`calculus::lunar::find_moon_altitude_range`] |
//! | [`Star`] | [`calculus::stellar::find_star_above_periods`] / [`calculus::stellar::find_star_range_periods`] |
//! | [`direction::ICRS`] | [`calculus::stellar::find_star_above_periods`] / [`calculus::stellar::find_star_range_periods`] |
//!
//! ## Quick Start
//!
//! ```rust
//! use siderust::bodies::Sun;
//! use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
//! use siderust::coordinates::centers::ObserverSite;
//! use siderust::coordinates::spherical::direction;
//! use siderust::time::{ModifiedJulianDate, Period};
//! use qtty::*;
//!
//! let site = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
//! let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
//!
//! let query = AltitudeQuery {
//!     observer: site,
//!     window,
//!     min_altitude: Degrees::new(-90.0),
//!     max_altitude: Degrees::new(0.0),
//! };
//!
//! // Same call shape for any body:
//! let sun_periods = Sun.altitude_periods(&query);
//! let star_periods = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716))
//!     .altitude_periods(&query);
//! ```

use super::types::AltitudeQuery;
use crate::bodies::solar_system;
use crate::bodies::Star;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::spherical::direction;
use crate::time::{ModifiedJulianDate, Period};
use qtty::*;

// ---------------------------------------------------------------------------
// Trait Definition
// ---------------------------------------------------------------------------

/// Unified interface for computing altitude periods of any celestial body.
///
/// Implementors delegate to the appropriate analytical/numerical engine in
/// the `calculus` layer.  The trait is intentionally small — one required
/// method plus convenience defaults.
pub trait AltitudePeriodsProvider {
    /// Returns all contiguous intervals inside `query.window` where the
    /// body's topocentric altitude is within
    /// `[query.min_altitude, query.max_altitude]`.
    ///
    /// The returned vector is sorted chronologically.  An empty vector
    /// means the body never enters the requested band during the window.
    fn altitude_periods(&self, query: &AltitudeQuery) -> Vec<Period<ModifiedJulianDate>>;

    /// Convenience: intervals where altitude is **above** `threshold`.
    ///
    /// Default implementation calls [`altitude_periods`](Self::altitude_periods)
    /// with `max_altitude = 90°`.
    fn above_threshold(
        &self,
        observer: ObserverSite,
        window: Period<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<Period<ModifiedJulianDate>> {
        self.altitude_periods(&AltitudeQuery {
            observer,
            window,
            min_altitude: threshold,
            max_altitude: Degrees::new(90.0),
        })
    }

    /// Convenience: intervals where altitude is **below** `threshold`.
    ///
    /// Default implementation calls [`altitude_periods`](Self::altitude_periods)
    /// with `min_altitude = −90°`.
    fn below_threshold(
        &self,
        observer: ObserverSite,
        window: Period<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<Period<ModifiedJulianDate>> {
        self.altitude_periods(&AltitudeQuery {
            observer,
            window,
            min_altitude: Degrees::new(-90.0),
            max_altitude: threshold,
        })
    }

    /// Compute the altitude of this body at a single instant.
    ///
    /// Returns the topocentric altitude in radians.
    fn altitude_at(&self, observer: &ObserverSite, mjd: ModifiedJulianDate) -> Radians;

    /// Hint for the scan step to use when searching for events.
    ///
    /// Returns `None` to use the default (10 minutes). Bodies with slower
    /// apparent motion (like the Moon) can return a larger step for efficiency.
    fn scan_step_hint(&self) -> Option<Days> {
        None
    }
}

// ---------------------------------------------------------------------------
// Free Function Entry Point
// ---------------------------------------------------------------------------

/// Generic entry point: compute altitude periods for any body that implements
/// [`AltitudePeriodsProvider`].
///
/// This is a thin wrapper around the trait method, provided for callers who
/// prefer a function‑style API.
///
/// # Example
/// ```rust
/// use siderust::calculus::altitude::{altitude_periods, AltitudeQuery};
/// use siderust::bodies::Sun;
/// use siderust::coordinates::centers::ObserverSite;
/// use siderust::time::{ModifiedJulianDate, Period};
/// use qtty::*;
///
/// let site = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
/// let query = AltitudeQuery {
///     observer: site,
///     window,
///     min_altitude: Degrees::new(0.0),
///     max_altitude: Degrees::new(90.0),
/// };
/// let periods = altitude_periods(&Sun, &query);
/// ```
#[inline]
pub fn altitude_periods<B: AltitudePeriodsProvider>(
    body: &B,
    query: &AltitudeQuery,
) -> Vec<Period<ModifiedJulianDate>> {
    body.altitude_periods(query)
}

// ---------------------------------------------------------------------------
// Implementations
// ---------------------------------------------------------------------------

/// **Sun** — delegates to [`calculus::solar`].
impl AltitudePeriodsProvider for solar_system::Sun {
    fn altitude_periods(&self, query: &AltitudeQuery) -> Vec<Period<ModifiedJulianDate>> {
        if query.window.duration_days() <= 0.0 {
            return Vec::new();
        }
        use crate::calculus::solar;

        // Fast path: full above query (max ≈ 90°)
        if query.max_altitude.value() >= 89.99 {
            solar::find_day_periods(query.observer, query.window, query.min_altitude)
        } else if query.min_altitude.value() <= -89.99 {
            // Full below query
            solar::find_night_periods(query.observer, query.window, query.max_altitude)
        } else {
            solar::find_sun_range_periods(
                query.observer,
                query.window,
                (query.min_altitude, query.max_altitude),
            )
        }
    }

    fn altitude_at(&self, observer: &ObserverSite, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::solar::sun_altitude_rad(mjd, observer)
    }
}

/// **Moon** — delegates to [`calculus::lunar`].
impl AltitudePeriodsProvider for solar_system::Moon {
    fn altitude_periods(&self, query: &AltitudeQuery) -> Vec<Period<ModifiedJulianDate>> {
        if query.window.duration_days() <= 0.0 {
            return Vec::new();
        }
        use crate::calculus::lunar;

        if query.max_altitude.value() >= 89.99 {
            lunar::find_moon_above_horizon(query.observer, query.window, query.min_altitude)
        } else if query.min_altitude.value() <= -89.99 {
            lunar::find_moon_below_horizon(query.observer, query.window, query.max_altitude)
        } else {
            lunar::find_moon_altitude_range(
                query.observer,
                query.window,
                (query.min_altitude, query.max_altitude),
            )
        }
    }

    fn altitude_at(&self, observer: &ObserverSite, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::lunar::moon_altitude_rad(mjd, observer)
    }

    fn scan_step_hint(&self) -> Option<Days> {
        // Moon moves slower, 2-hour steps are sufficient
        Some(Hours::new(2.0).to::<Day>())
    }
}

/// **Star** — extracts RA/Dec from the star's target, delegates to
/// [`calculus::stellar`].
impl AltitudePeriodsProvider for Star<'_> {
    fn altitude_periods(&self, query: &AltitudeQuery) -> Vec<Period<ModifiedJulianDate>> {
        let dir = direction::ICRS::from(self);
        dir.altitude_periods(query)
    }

    fn altitude_at(&self, observer: &ObserverSite, mjd: ModifiedJulianDate) -> Radians {
        let dir = direction::ICRS::from(self);
        dir.altitude_at(observer, mjd)
    }
}

/// **direction::ICRS** — the lightest path: raw RA/Dec → stellar engine.
impl AltitudePeriodsProvider for direction::ICRS {
    fn altitude_periods(&self, query: &AltitudeQuery) -> Vec<Period<ModifiedJulianDate>> {
        if query.window.duration_days() <= 0.0 {
            return Vec::new();
        }
        use crate::calculus::stellar;

        if query.max_altitude.value() >= 89.99 {
            stellar::find_star_above_periods(
                self.ra(),
                self.dec(),
                query.observer,
                query.window,
                query.min_altitude,
            )
        } else if query.min_altitude.value() <= -89.99 {
            stellar::find_star_below_periods(
                self.ra(),
                self.dec(),
                query.observer,
                query.window,
                query.max_altitude,
            )
        } else {
            stellar::find_star_range_periods(
                self.ra(),
                self.dec(),
                query.observer,
                query.window,
                (query.min_altitude, query.max_altitude),
            )
        }
    }

    fn altitude_at(&self, observer: &ObserverSite, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::stellar::fixed_star_altitude_rad(mjd, observer, self.ra(), self.dec())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog;

    fn greenwich() -> ObserverSite {
        ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn one_day_window() -> Period<ModifiedJulianDate> {
        Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60001.0),
        )
    }

    fn one_week_window() -> Period<ModifiedJulianDate> {
        Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60007.0),
        )
    }

    // --- Consistent API shape ---

    #[test]
    fn sun_above_horizon_via_trait() {
        let periods =
            solar_system::Sun.above_threshold(greenwich(), one_day_window(), Degrees::new(0.0));
        assert!(!periods.is_empty(), "Sun should be above horizon at 51°N");
        for p in &periods {
            assert!(p.duration_days() > 0.0);
            assert!(p.duration_days() < 1.0);
        }
    }

    #[test]
    fn moon_above_horizon_via_trait() {
        let periods =
            solar_system::Moon.above_threshold(greenwich(), one_week_window(), Degrees::new(0.0));
        assert!(
            !periods.is_empty(),
            "Moon should be above horizon at some point in a week"
        );
    }

    #[test]
    fn star_above_horizon_via_trait() {
        let sirius = &catalog::SIRIUS;
        let periods = sirius.above_threshold(greenwich(), one_day_window(), Degrees::new(0.0));
        // Sirius (Dec ≈ −16.7°) rises and sets at 51°N
        assert!(
            !periods.is_empty(),
            "Sirius should be above horizon for part of the day"
        );
    }

    #[test]
    fn icrs_direction_above_horizon_via_trait() {
        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let periods = sirius_dir.above_threshold(greenwich(), one_day_window(), Degrees::new(0.0));
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

        let star_periods = sirius.above_threshold(observer, window, Degrees::new(0.0));
        let dir_periods = sirius_dir.above_threshold(observer, window, Degrees::new(0.0));

        assert_eq!(
            star_periods.len(),
            dir_periods.len(),
            "Star and direction::ICRS should produce the same number of periods"
        );
        for (sp, dp) in star_periods.iter().zip(dir_periods.iter()) {
            assert!(
                (sp.start.value() - dp.start.value()).abs() < 1e-6,
                "Period starts should match"
            );
            assert!(
                (sp.end.value() - dp.end.value()).abs() < 1e-6,
                "Period ends should match"
            );
        }
    }

    // --- Generic free function ---

    #[test]
    fn free_function_works_for_sun() {
        let query = AltitudeQuery {
            observer: greenwich(),
            window: one_day_window(),
            min_altitude: Degrees::new(0.0),
            max_altitude: Degrees::new(90.0),
        };
        let periods = altitude_periods(&solar_system::Sun, &query);
        assert!(!periods.is_empty());
    }

    #[test]
    fn free_function_works_for_icrs_direction() {
        let dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let query = AltitudeQuery {
            observer: greenwich(),
            window: one_day_window(),
            min_altitude: Degrees::new(0.0),
            max_altitude: Degrees::new(90.0),
        };
        let periods = altitude_periods(&dir, &query);
        assert!(!periods.is_empty());
    }

    // --- altitude_at single-point ---

    #[test]
    fn altitude_at_consistent_across_types() {
        let observer = greenwich();
        let mjd = ModifiedJulianDate::new(51544.5); // J2000 epoch in MJD

        let sun_alt = solar_system::Sun.altitude_at(&observer, mjd);
        assert!(sun_alt.value().abs() < std::f64::consts::FRAC_PI_2);

        let moon_alt = solar_system::Moon.altitude_at(&observer, mjd);
        assert!(moon_alt.value().abs() < std::f64::consts::FRAC_PI_2);

        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let star_alt = sirius_dir.altitude_at(&observer, mjd);
        assert!(star_alt.value().abs() < std::f64::consts::FRAC_PI_2);
    }

    // --- Edge cases ---

    #[test]
    fn full_sky_range_returns_full_window() {
        let query = AltitudeQuery {
            observer: greenwich(),
            window: one_day_window(),
            min_altitude: Degrees::new(-90.0),
            max_altitude: Degrees::new(90.0),
        };
        let periods = solar_system::Sun.altitude_periods(&query);
        // The sun's altitude is always in [-90, 90], so we should get the full window
        assert!(
            !periods.is_empty(),
            "Full sky range should return at least one period"
        );
        let total: f64 = periods.iter().map(|p| p.duration_days()).sum();
        assert!(
            (total - 1.0).abs() < 0.01,
            "Full sky range should span ~1 day, got {} days",
            total
        );
    }

    #[test]
    fn polaris_circumpolar_via_trait() {
        let polaris = &catalog::POLARIS;
        let periods = polaris.above_threshold(greenwich(), one_day_window(), Degrees::new(0.0));
        assert_eq!(
            periods.len(),
            1,
            "Polaris should be continuously above horizon at 51°N"
        );
        assert!(
            (periods[0].duration_days() - 1.0).abs() < 0.01,
            "Polaris up-period should span the full day"
        );
    }

    #[test]
    fn polaris_never_below_minus90_via_trait() {
        let polaris = &catalog::POLARIS;
        // Polaris is circumpolar at 51°N — should never be below -90° (vacuous)
        let periods = polaris.below_threshold(greenwich(), one_day_window(), Degrees::new(-80.0));
        assert!(
            periods.is_empty(),
            "Polaris should never be below -80° at 51°N"
        );
    }

    #[test]
    fn empty_window_returns_empty() {
        let window = Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60000.0),
        );
        let query = AltitudeQuery {
            observer: greenwich(),
            window,
            min_altitude: Degrees::new(0.0),
            max_altitude: Degrees::new(90.0),
        };
        let periods = solar_system::Sun.altitude_periods(&query);
        assert!(periods.is_empty(), "Empty window should return no periods");
    }

    #[test]
    fn below_threshold_sun_night_via_trait() {
        let nights =
            solar_system::Sun.below_threshold(greenwich(), one_week_window(), Degrees::new(-18.0));
        assert!(!nights.is_empty(), "Should find astronomical night at 51°N");
    }

    #[test]
    fn altitude_range_twilight_via_trait() {
        let query = AltitudeQuery {
            observer: greenwich(),
            window: Period::new(
                ModifiedJulianDate::new(60000.0),
                ModifiedJulianDate::new(60002.0),
            ),
            min_altitude: Degrees::new(-18.0),
            max_altitude: Degrees::new(-12.0),
        };
        let bands = solar_system::Sun.altitude_periods(&query);
        assert!(
            bands.len() >= 2,
            "Should find at least 2 twilight bands in 2 days, found {}",
            bands.len()
        );
    }

    // --- Periods are sorted and non-overlapping ---

    #[test]
    fn periods_are_sorted_and_non_overlapping() {
        let sirius_dir = direction::ICRS::new(Degrees::new(101.287), Degrees::new(-16.716));
        let periods = sirius_dir.above_threshold(greenwich(), one_week_window(), Degrees::new(0.0));
        for w in periods.windows(2) {
            assert!(
                w[0].end.value() <= w[1].start.value(),
                "Periods should be non-overlapping and sorted: {:?} vs {:?}",
                w[0],
                w[1]
            );
        }
    }
}
