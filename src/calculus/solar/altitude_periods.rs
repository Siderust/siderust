//! # Solar Altitude Window Periods
//!
//! Sun-specific convenience wrappers around the generic altitude-window routines
//! in [`crate::calculus::events::altitude_periods`].
//!
//! This module exists to keep `calculus::events` generic (body-agnostic) while
//! still providing ergonomic helpers for common solar concepts like twilight.

use crate::astro::JulianDate;
use crate::bodies::solar_system::Sun;
use crate::calculus::events::altitude_periods::{find_altitude_periods, AltitudeCondition, AltitudePeriod};
use crate::coordinates::centers::ObserverSite;
use crate::time::{ModifiedJulianDate, Period};
use qtty::{AstronomicalUnit, Degrees, Radian};

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
/// Positive above the horizon, negative below.
pub fn sun_altitude_rad(jd: JulianDate, site: &ObserverSite) -> f64 {
    let horiz = Sun::get_horizontal::<AstronomicalUnit>(jd, *site);
    horiz.alt().to::<Radian>().value()
}

/// Common twilight types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Twilight {
    Civil,
    Nautical,
    Astronomical,
    Horizon,
    ApparentHorizon,
}

impl From<Twilight> for Degrees {
    fn from(t: Twilight) -> Degrees {
        match t {
            Twilight::Civil => Degrees::new(-6.0),
            Twilight::Nautical => Degrees::new(-12.0),
            Twilight::Astronomical => Degrees::new(-18.0),
            Twilight::Horizon => Degrees::new(0.0),
            Twilight::ApparentHorizon => Degrees::new(-0.833),
        }
    }
}

/// Finds night periods (Sun below `twilight`) inside `period`.
pub fn find_night_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    let tw: Degrees = twilight.into();
    find_altitude_periods(altitude_fn, period, AltitudeCondition::below(tw))
}

/// Finds day periods (Sun above `twilight`) inside `period`.
pub fn find_day_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    let tw: Degrees = twilight.into();
    find_altitude_periods(altitude_fn, period, AltitudeCondition::above(tw))
}

/// Backwards-compatible alias expected by `calculus::solar::mod.rs`.
pub fn find_sun_above_altitude<T: Into<Degrees>>(site: ObserverSite, period: Period<ModifiedJulianDate>, twilight: T) -> Option<Vec<AltitudePeriod>> {
    find_day_periods(site, period, twilight)
}

/// Finds periods where Sun altitude is within `range` (min, max) inside `period`.
pub fn find_sun_range_periods(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    find_altitude_periods(altitude_fn, period, AltitudeCondition::between(range.0, range.1))
}

/// Backwards-compatible alias expected by `calculus::solar::mod.rs`.
pub fn find_sun_in_altitude_range(site: ObserverSite, period: Period<ModifiedJulianDate>, range: (Degrees, Degrees)) -> Option<Vec<AltitudePeriod>> {
    find_sun_range_periods(site, period, range)
}

/// Standard twilight threshold definitions (Sun center altitude).
pub mod twilight {
    use qtty::Degrees;

    /// Civil twilight: Sun center 6° below horizon (-6°)
    pub const CIVIL: Degrees = Degrees::new(-6.0);

    /// Nautical twilight: Sun center 12° below horizon (-12°)
    pub const NAUTICAL: Degrees = Degrees::new(-12.0);

    /// Astronomical twilight: Sun center 18° below horizon (-18°)
    pub const ASTRONOMICAL: Degrees = Degrees::new(-18.0);

    /// Sunrise/sunset: Sun center at geometric horizon (0°)
    /// Note: For apparent sunrise/sunset, use -0.833° to account for refraction
    pub const HORIZON: Degrees = Degrees::new(0.0);

    /// Apparent sunrise/sunset accounting for atmospheric refraction (-0.833°)
    pub const APPARENT_HORIZON: Degrees = Degrees::new(-0.833);
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    fn greenwich_site() -> ObserverSite {
        ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Quantity::<Meter>::new(0.0))
    }

    #[test]
    fn test_sun_altitude_basic() {
        let site = greenwich_site();
        let jd = JulianDate::J2000;
        let alt = sun_altitude_rad(jd, &site);
        assert!(alt > -std::f64::consts::FRAC_PI_2 && alt < std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_find_night_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let period = Period::new(mjd_start, mjd_end);

        let nights = find_night_periods(site, period, twilight::ASTRONOMICAL);
        assert!(nights.is_some(), "Should find night periods at 51° latitude");

        let nights = nights.unwrap();
        assert!(!nights.is_empty(), "Should have at least one night period");

        for night in &nights {
            assert!(night.duration_days() > 0.0, "Night duration should be positive");
            assert!(night.duration_days() < 1.0, "Night should be less than 24 hours");
        }
    }

    #[test]
    fn test_find_altitude_range_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);

        let period = Period::new(mjd_start, mjd_end);

        let nights = find_sun_range_periods(
            site,
            period,
            (Degrees::new(-90.0), Degrees::new(-18.0)),
        );

        assert!(nights.is_some(), "Should find night periods using range");
        let nights = nights.unwrap();
        assert!(!nights.is_empty(), "Should have at least one night period");

        let nautical = find_sun_range_periods(
            site,
            period,
            (Degrees::new(-18.0), Degrees::new(-12.0)),
        );

        assert!(nautical.is_some(), "Should find nautical twilight periods");
    }
}
