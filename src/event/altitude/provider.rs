// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Provider, Trait-Based Dispatch Layer
//!
//! Defines [`AltitudeProvider`] and implementations for celestial targets.

use crate::bodies::solar_system;
use crate::bodies::Star;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::spherical::direction;
use crate::qtty::*;
use crate::time::ModifiedJulianDate;

use crate::coordinates::{cartesian, centers::Geocentric, frames};
use crate::event::horizontal;
use crate::time::JulianDate;

/// Unified interface for evaluating topocentric altitude of any celestial target.
///
/// Implementors delegate single-point altitude to the appropriate analytical or
/// numerical engine. Period queries use [`super::above_threshold`],
/// [`super::below_threshold`], and [`super::altitude_ranges`].
pub trait AltitudeProvider {
    /// Compute the altitude of this body at a single instant (radians).
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
    fn scan_step_hint(&self) -> Option<Days> {
        None
    }
}

impl AltitudeProvider for solar_system::Sun {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::event::solar::sun_altitude_rad(mjd, observer)
    }
}

impl AltitudeProvider for solar_system::Moon {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        crate::event::lunar::moon_altitude_rad(mjd, observer)
    }

    fn scan_step_hint(&self) -> Option<Days> {
        Some(Hours::new(2.0).to::<Day>())
    }
}

impl AltitudeProvider for Star<'_> {
    fn altitude_at(&self, observer: &Geodetic<ECEF>, mjd: ModifiedJulianDate) -> Radians {
        direction::ICRS::from(self).altitude_at(observer, mjd)
    }

    fn altitude_at_with_policy(
        &self,
        observer: &Geodetic<ECEF>,
        mjd: ModifiedJulianDate,
        policy: crate::astro::apparent::CorrectionPolicy,
    ) -> Radians {
        direction::ICRS::from(self).altitude_at_with_policy(observer, mjd, policy)
    }
}

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
}

use crate::coordinates::transform::Transform;

const PLANET_SCAN_STEP: Days = Quantity::<Hour>::new(2.0).to_const::<Day>();

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
    let bary_ecl = vsop87e_fn(jd);
    let geo_equ: cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, AstronomicalUnit> =
        bary_ecl.transform(jd);
    let topo = horizontal::geocentric_j2000_to_apparent_topocentric(&geo_equ, *site, jd);
    let horiz = horizontal::equatorial_to_horizontal(&topo, *site, jd);
    horiz.alt().to::<Radian>()
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog;
    use crate::event::altitude::{above_threshold, altitude_ranges, below_threshold, SearchOpts};
    use crate::time::Interval;

    struct FixedAltitudeTarget {
        alt: Radians,
    }

    impl AltitudeProvider for FixedAltitudeTarget {
        fn altitude_at(&self, _observer: &Geodetic<ECEF>, _mjd: ModifiedJulianDate) -> Radians {
            self.alt
        }
    }

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
    fn custom_altitude_provider_uses_generic_search_path() {
        let target = FixedAltitudeTarget {
            alt: Degrees::new(45.0).to::<Radian>(),
        };
        let periods = above_threshold(
            &target,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert_eq!(periods.len(), 1);
        assert!(
            ((periods[0].end.raw() - periods[0].start.raw()) - Days::new(1.0)).abs()
                < Days::new(0.01)
        );
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
        assert!(!periods.is_empty());
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
        assert!(!periods.is_empty());
    }

    #[test]
    fn star_above_horizon_via_standard_api() {
        let periods = above_threshold(
            &catalog::SIRIUS,
            &greenwich(),
            one_day_window(),
            Degrees::new(0.0),
            SearchOpts::default(),
        );
        assert!(!periods.is_empty());
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
        assert!(!nights.is_empty());
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
        assert!(bands.len() >= 2);
    }
}
