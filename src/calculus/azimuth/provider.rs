// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Azimuth Provider — Trait-Based Dispatch Layer
//!
//! Defines [`AzimuthProvider`] and implementations that normalise the azimuth
//! API across Sun, Moon, and fixed stars.
//!
//! ## Design
//!
//! The trait is the *dispatch* layer — all astronomical math lives in the
//! per-body engine modules (`calculus::solar`, `calculus::lunar`,
//! `calculus::stellar`).  Each `impl` simply delegates.
//!
//! | Body type | Engine |
//! |-----------|--------|
//! | [`solar_system::Sun`]  | [`calculus::solar::sun_azimuth_rad`]         |
//! | [`solar_system::Moon`] | [`calculus::lunar::moon_azimuth_rad`]        |
//! | [`Star`]               | delegates to [`direction::ICRS`] impl below  |
//! | [`direction::ICRS`]    | [`calculus::stellar::fixed_star_azimuth_rad`] |
//!
//! ## Quick Start
//!
//! ```rust
//! use siderust::bodies::Sun;
//! use siderust::calculus::azimuth::{AzimuthProvider, AzimuthQuery};
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::ECEF;
//! use siderust::time::{ModifiedJulianDate, MJD, Period};
//! use qtty::*;
//!
//! let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
//! let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
//! let query = AzimuthQuery {
//!     observer: site,
//!     window,
//!     min_azimuth: Degrees::new(90.0),
//!     max_azimuth: Degrees::new(270.0),
//! };
//!
//! // Same call shape for any body:
//! let sun_periods = Sun.azimuth_periods(&query);
//! ```

use super::events;
use super::types::AzimuthQuery;
use crate::bodies::solar_system;
use crate::bodies::Star;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::spherical::direction;
use crate::time::{complement_within, ModifiedJulianDate, Period, MJD};
use qtty::*;

// ---------------------------------------------------------------------------
// Trait Definition
// ---------------------------------------------------------------------------

/// Unified interface for computing azimuth-related quantities and periods
/// for any celestial body.
///
/// ## Time Scale
///
/// All `ModifiedJulianDate` and `Period<MJD>` values are on the canonical
/// JD(TT) axis.  Convert UTC instants with `ModifiedJulianDate::from_utc(…)`
/// before using this API.
///
/// ## Azimuth convention
///
/// North-clockwise (North = 0°), values in `[0°, 360°)`.
pub trait AzimuthProvider {
    /// Compute the topocentric azimuth of this body at a single instant.
    ///
    /// Returns a value in radians, range `[0, 2π)` (North-clockwise).
    fn azimuth_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians;

    /// Returns all contiguous intervals inside `query.window` where the body's
    /// azimuth is within `[query.min_azimuth, query.max_azimuth]`.
    ///
    /// If `min_azimuth > max_azimuth` the query is a **wrap-around** (North-
    /// crossing) range and the result covers the combined arc.
    ///
    /// The returned vector is sorted chronologically.  An empty vector means
    /// the body never enters the requested band during the window.
    fn azimuth_periods(&self, query: &AzimuthQuery) -> Vec<Period<MJD>>;

    /// Convenience: intervals where azimuth is within `[min_az, max_az]`.
    ///
    /// Equivalent to calling [`azimuth_periods`](Self::azimuth_periods)
    /// directly.
    fn in_azimuth_range(
        &self,
        observer: Geodetic<ECEF>,
        window: Period<MJD>,
        min_az: Degrees,
        max_az: Degrees,
    ) -> Vec<Period<MJD>> {
        self.azimuth_periods(&AzimuthQuery {
            observer,
            window,
            min_azimuth: min_az,
            max_azimuth: max_az,
        })
    }

    /// Convenience: intervals where azimuth is **outside** `[min_az, max_az]`.
    ///
    /// Returns the complement of [`in_azimuth_range`](Self::in_azimuth_range)
    /// within `window`.
    fn outside_azimuth_range(
        &self,
        observer: Geodetic<ECEF>,
        window: Period<MJD>,
        min_az: Degrees,
        max_az: Degrees,
    ) -> Vec<Period<MJD>> {
        let inside = self.in_azimuth_range(observer, window, min_az, max_az);
        complement_within(window, &inside)
    }

    /// Hint for the scan step to use when searching for azimuth events.
    ///
    /// Returns `None` to use the default (10 minutes). Bodies with slower
    /// apparent motion (like the Moon) can return a larger step.
    fn scan_step_hint(&self) -> Option<Days> {
        None
    }
}

// ---------------------------------------------------------------------------
// Free Function Entry Point
// ---------------------------------------------------------------------------

/// Generic entry point: compute azimuth periods for any body implementing
/// [`AzimuthProvider`].
///
/// This is a thin wrapper around the trait method for callers preferring a
/// function-style API.
///
/// # Example
/// ```rust
/// use siderust::calculus::azimuth::{azimuth_periods, AzimuthQuery};
/// use siderust::bodies::Sun;
/// use siderust::coordinates::centers::Geodetic;
/// use siderust::coordinates::frames::ECEF;
/// use siderust::time::{ModifiedJulianDate, MJD, Period};
/// use qtty::*;
///
/// let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let window = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60001.0));
/// let query = AzimuthQuery {
///     observer: site,
///     window,
///     min_azimuth: Degrees::new(90.0),
///     max_azimuth: Degrees::new(270.0),
/// };
/// let periods = azimuth_periods(&Sun, &query);
/// ```
#[inline]
pub fn azimuth_periods<B: AzimuthProvider>(body: &B, query: &AzimuthQuery) -> Vec<Period<MJD>> {
    body.azimuth_periods(query)
}

// ---------------------------------------------------------------------------
// Implementations
// ---------------------------------------------------------------------------

/// **Sun** — delegates to [`calculus::solar::sun_azimuth_rad`].
impl AzimuthProvider for solar_system::Sun {
    fn azimuth_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::solar::sun_azimuth_rad(mjd, observer)
    }

    fn azimuth_periods(&self, query: &AzimuthQuery) -> Vec<Period<MJD>> {
        if query.window.duration() <= Days::zero() {
            return Vec::new();
        }
        events::azimuth_range_periods(self, query)
    }

    fn scan_step_hint(&self) -> Option<Days> {
        Some(Hours::new(2.0).to::<Day>())
    }
}

/// **Moon** — delegates to [`calculus::lunar::moon_azimuth_rad`].
impl AzimuthProvider for solar_system::Moon {
    fn azimuth_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::lunar::moon_azimuth_rad(mjd, observer)
    }

    fn azimuth_periods(&self, query: &AzimuthQuery) -> Vec<Period<MJD>> {
        if query.window.duration() <= Days::zero() {
            return Vec::new();
        }
        events::azimuth_range_periods(self, query)
    }

    fn scan_step_hint(&self) -> Option<Days> {
        // Moon moves slower; 2-hour steps are sufficient.
        Some(Hours::new(2.0).to::<Day>())
    }
}

/// **Star** — extracts RA/Dec, delegates to [`direction::ICRS`].
impl AzimuthProvider for Star<'_> {
    fn azimuth_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        let dir = direction::ICRS::from(self);
        dir.azimuth_at(observer, mjd)
    }

    fn azimuth_periods(&self, query: &AzimuthQuery) -> Vec<Period<MJD>> {
        let dir = direction::ICRS::from(self);
        dir.azimuth_periods(query)
    }
}

/// **direction::ICRS** — raw RA/Dec → stellar engine.
impl AzimuthProvider for direction::ICRS {
    fn azimuth_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::calculus::stellar::fixed_star_azimuth_rad(mjd, observer, self.ra(), self.dec())
    }

    fn azimuth_periods(&self, query: &AzimuthQuery) -> Vec<Period<MJD>> {
        if query.window.duration() <= Days::zero() {
            return Vec::new();
        }
        events::azimuth_range_periods(self, query)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog;
    use crate::bodies::solar_system::Moon;

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    fn one_day_window() -> Period<MJD> {
        Period::new(
            ModifiedJulianDate::new(60000.0),
            ModifiedJulianDate::new(60001.0),
        )
    }

    #[test]
    fn sun_azimuth_at_returns_valid_range() {
        let az = solar_system::Sun.azimuth_at(
            &greenwich(),
            ModifiedJulianDate::new(60000.5), // noon-ish
        );
        assert!(az.value() >= 0.0, "azimuth must be ≥ 0");
        assert!(az.value() < std::f64::consts::TAU, "azimuth must be < 2π");
    }

    #[test]
    fn moon_azimuth_at_returns_valid_range() {
        let az = Moon.azimuth_at(&greenwich(), ModifiedJulianDate::new(60000.5));
        assert!(az.value() >= 0.0);
        assert!(az.value() < std::f64::consts::TAU);
    }

    #[test]
    fn star_azimuth_at_returns_valid_range() {
        let sirius = &catalog::SIRIUS;
        let az = sirius.azimuth_at(&greenwich(), ModifiedJulianDate::new(60000.5));
        assert!(az.value() >= 0.0);
        assert!(az.value() < std::f64::consts::TAU);
    }

    #[test]
    fn star_and_icrs_agree() {
        let sirius = &catalog::SIRIUS;
        let dir = direction::ICRS::from(sirius);
        let mjd = ModifiedJulianDate::new(60000.5);
        let az_star = sirius.azimuth_at(&greenwich(), mjd);
        let az_dir = dir.azimuth_at(&greenwich(), mjd);
        assert!(
            (az_star.value() - az_dir.value()).abs() < 1e-10,
            "Star and direction::ICRS must agree"
        );
    }

    #[test]
    fn sun_azimuth_periods_eastern_half() {
        let query = AzimuthQuery {
            observer: greenwich(),
            window: one_day_window(),
            min_azimuth: Degrees::new(90.0),
            max_azimuth: Degrees::new(270.0),
        };
        let periods = solar_system::Sun.azimuth_periods(&query);
        // Sun is on the eastern half (az 90-270) for roughly half a day
        assert!(
            !periods.is_empty(),
            "Sun should cross the eastern half in 24h"
        );
    }

    #[test]
    fn outside_azimuth_range_is_complement() {
        let observer = greenwich();
        let window = one_day_window();
        let inside = solar_system::Sun.in_azimuth_range(
            observer,
            window,
            Degrees::new(90.0),
            Degrees::new(270.0),
        );
        let outside = solar_system::Sun.outside_azimuth_range(
            observer,
            window,
            Degrees::new(90.0),
            Degrees::new(270.0),
        );

        // Total covered should be the full window
        let total_inside: f64 = inside.iter().map(|p| p.duration_days().value()).sum();
        let total_outside: f64 = outside.iter().map(|p| p.duration_days().value()).sum();
        let window_len = window.duration_days().value();
        assert!(
            (total_inside + total_outside - window_len).abs() < 1e-6,
            "inside + outside should equal the full window"
        );
    }
}
