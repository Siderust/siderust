// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::Moon;

use crate::calculus::ephemeris::Ephemeris;
use crate::calculus::horizontal;
use crate::calculus::lunar::phase::{
    find_phase_events, illumination_above, illumination_below, illumination_range,
    moon_phase_geocentric, moon_phase_topocentric, MoonPhaseGeometry, PhaseEvent, PhaseSearchOpts,
};
use crate::coordinates::transform::context::DefaultEphemeris;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use crate::time::{JulianDate, Period, MJD};
use qtty::{AstronomicalUnits, Kilometer, LengthUnit, Meter, Quantity};

impl Moon {
    /// Returns the **apparent topocentric equatorial coordinates** of the Moon
    /// as seen from a given `Geodetic<ECEF>` at the specified Julian Date.
    ///
    /// This method accounts for:
    /// - **Topocentric parallax**: Critical for the Moon due to its proximity (~1° at horizon)
    /// - **Precession**: J2000 → mean-of-date transformation
    /// - **Nutation**: Mean-of-date → true-of-date correction
    ///
    /// ### Parameters
    /// - `jd`: Julian Day for which to compute the Moon's apparent position
    /// - `site`: Observer location on Earth
    ///
    /// ### Returns
    /// A `spherical::Position<Topocentric, EquatorialTrueOfDate, U>` representing the Moon's
    /// apparent right ascension and declination from the observer's location.
    ///
    /// ### Notes
    /// Unlike the Sun, topocentric parallax correction is **essential** for the Moon
    /// due to its proximity to Earth (average distance ~384,400 km).
    pub fn get_apparent_topocentric_equ<U: LengthUnit>(
        jd: JulianDate,
        site: Geodetic<frames::ECEF>,
    ) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<Quantity<Kilometer>> + From<AstronomicalUnits>,
    {
        // 1) Get Moon's geocentric ecliptic position from the active ephemeris backend
        let moon_geo_ecliptic: cartesian::Position<
            Geocentric,
            frames::EclipticMeanJ2000,
            Kilometer,
        > = DefaultEphemeris::moon_geocentric(jd);

        // 2) Transform: EclipticMeanJ2000 → EquatorialMeanJ2000
        let moon_geo_eq_j2000: cartesian::Position<
            Geocentric,
            frames::EquatorialMeanJ2000,
            Kilometer,
        > = TransformFrame::to_frame(&moon_geo_ecliptic);

        // 3-5) Shared pipeline in Kilometer: topocentric parallax → precession → nutation → spherical
        let topo_sph_km: spherical::Position<Topocentric, frames::EquatorialTrueOfDate, Kilometer> =
            horizontal::geocentric_j2000_to_apparent_topocentric::<Kilometer>(
                &moon_geo_eq_j2000,
                site,
                jd,
            );

        // 6) Convert from Kilometer to target unit U
        let dist_u: Quantity<U> = topo_sph_km.distance.into();

        affn::spherical::Position::<Topocentric, frames::EquatorialTrueOfDate, U>::new_raw_with_params(
            *topo_sph_km.center_params(),
            topo_sph_km.polar,
            topo_sph_km.azimuth,
            dist_u,
        )
    }

    /// Returns the Moon's **horizontal coordinates** (altitude, azimuth) as seen
    /// from a given `Geodetic<ECEF>` at the specified Julian Date or Modified Julian Date.
    ///
    /// This is a convenience wrapper that computes the apparent topocentric equatorial
    /// position and transforms it to horizontal coordinates.
    ///
    /// ### Parameters
    /// - `time`: Any type that can be converted to `JulianDate` (JD or Mjd)
    /// - `site`: Observer location on Earth
    ///
    /// ### Returns
    /// A `spherical::Position<Topocentric, Horizontal, U>` with:
    /// - Altitude (polar): elevation above horizon in degrees, [-90°, +90°]
    /// - Azimuth: bearing from North through East in degrees, [0°, 360°)
    /// - Distance: in the specified length unit
    ///
    /// ### Example
    /// ```rust
    /// use siderust::bodies::solar_system::Moon;
    /// use siderust::coordinates::centers::Geodetic;
    /// use siderust::coordinates::frames::ECEF;
    /// use siderust::time::{JulianDate, ModifiedJulianDate};
    /// use qtty::*;
    ///
    /// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Using JulianDate
    /// let moon_pos = Moon::get_horizontal::<Kilometer>(JulianDate::J2000, site);
    /// println!("Moon altitude: {}", moon_pos.alt().to::<Deg>());
    ///
    /// // Using ModifiedJulianDate
    /// let mjd = ModifiedJulianDate::new(60000.0);
    /// let moon_pos = Moon::get_horizontal::<Kilometer>(mjd, site);
    /// ```
    pub fn get_horizontal<U: LengthUnit>(
        time: impl Into<JulianDate>,
        site: Geodetic<frames::ECEF>,
    ) -> spherical::Position<Topocentric, frames::Horizontal, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<Quantity<Kilometer>> + From<AstronomicalUnits>,
    {
        let jd = time.into();
        let eq = Self::get_apparent_topocentric_equ::<U>(jd, site);
        horizontal::equatorial_to_horizontal(&eq, site, jd)
    }
}

// ===========================================================================
// Moon — phase API (DefaultEphemeris convenience methods)
// ===========================================================================

impl Moon {
    /// Geocentric Moon phase geometry at `time`, using the compile-time
    /// `DefaultEphemeris` backend.
    ///
    /// Returns illuminated fraction, phase angle, elongation, and waxing flag.
    /// For a specific ephemeris backend use
    /// [`moon_phase_geocentric::<E>(jd)`](crate::calculus::lunar::phase::moon_phase_geocentric).
    ///
    /// # Example
    /// ```rust
    /// use siderust::bodies::solar_system::Moon;
    /// use siderust::time::JulianDate;
    ///
    /// let geom = Moon::phase_geocentric(JulianDate::J2000);
    /// assert!(geom.illuminated_fraction >= 0.0 && geom.illuminated_fraction <= 1.0);
    /// ```
    pub fn phase_geocentric(time: impl Into<JulianDate>) -> MoonPhaseGeometry {
        moon_phase_geocentric::<DefaultEphemeris>(time.into())
    }

    /// Topocentric Moon phase geometry at `time` for the given observer `site`,
    /// using the compile-time `DefaultEphemeris` backend.
    ///
    /// For a specific ephemeris backend use
    /// [`moon_phase_topocentric::<E>(jd, site)`](crate::calculus::lunar::phase::moon_phase_topocentric).
    ///
    /// # Example
    /// ```rust
    /// use siderust::bodies::solar_system::Moon;
    /// use siderust::coordinates::centers::Geodetic;
    /// use siderust::coordinates::frames::ECEF;
    /// use siderust::time::JulianDate;
    /// use qtty::*;
    ///
    /// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.48 * DEG, 0.0 * M);
    /// let geom = Moon::phase_topocentric(JulianDate::J2000, site);
    /// println!("Illuminated: {:.1} %", geom.illuminated_percent());
    /// ```
    pub fn phase_topocentric(
        time: impl Into<JulianDate>,
        site: Geodetic<frames::ECEF>,
    ) -> MoonPhaseGeometry {
        moon_phase_topocentric::<DefaultEphemeris>(time.into(), site)
    }

    /// Find principal phase events (New Moon, First Quarter, Full Moon, Last
    /// Quarter) inside `window`, using the compile-time `DefaultEphemeris`.
    ///
    /// For a specific ephemeris backend use
    /// [`find_phase_events::<E>(window, opts)`](crate::calculus::lunar::phase::find_phase_events).
    ///
    /// # Example
    /// ```rust
    /// use siderust::bodies::solar_system::Moon;
    /// use siderust::calculus::lunar::phase::PhaseSearchOpts;
    /// use siderust::time::{JulianDate, ModifiedJulianDate, Period};
    /// use qtty::Days;
    ///
    /// let start  = ModifiedJulianDate::from(JulianDate::J2000);
    /// let window = Period::new(start, start + Days::new(35.0));
    /// let events = Moon::phase_events(window, PhaseSearchOpts::default());
    /// assert!(!events.is_empty());
    /// ```
    pub fn phase_events(window: Period<MJD>, opts: PhaseSearchOpts) -> Vec<PhaseEvent> {
        find_phase_events::<DefaultEphemeris>(window, opts)
    }

    /// Time windows inside `window` where geocentric illuminated fraction ≥
    /// `k_min`, using the compile-time `DefaultEphemeris`.
    pub fn illumination_above(
        window: Period<MJD>,
        k_min: f64,
        opts: PhaseSearchOpts,
    ) -> Vec<Period<MJD>> {
        illumination_above::<DefaultEphemeris>(window, k_min, opts)
    }

    /// Time windows inside `window` where geocentric illuminated fraction ≤
    /// `k_max`, using the compile-time `DefaultEphemeris`.
    pub fn illumination_below(
        window: Period<MJD>,
        k_max: f64,
        opts: PhaseSearchOpts,
    ) -> Vec<Period<MJD>> {
        illumination_below::<DefaultEphemeris>(window, k_max, opts)
    }

    /// Time windows inside `window` where geocentric illuminated fraction is
    /// within `[k_min, k_max]`, using the compile-time `DefaultEphemeris`.
    ///
    /// # Example
    /// ```rust
    /// use siderust::bodies::solar_system::Moon;
    /// use siderust::calculus::lunar::phase::PhaseSearchOpts;
    /// use siderust::time::{JulianDate, ModifiedJulianDate, Period};
    /// use qtty::Days;
    ///
    /// let start  = ModifiedJulianDate::from(JulianDate::J2000);
    /// let window = Period::new(start, start + Days::new(30.0));
    /// // Crescent phase: 5–35% illuminated
    /// let crescent = Moon::illumination_range(window, 0.05, 0.35, PhaseSearchOpts::default());
    /// ```
    pub fn illumination_range(
        window: Period<MJD>,
        k_min: f64,
        k_max: f64,
        opts: PhaseSearchOpts,
    ) -> Vec<Period<MJD>> {
        illumination_range::<DefaultEphemeris>(window, k_min, k_max, opts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::Geodetic;
    use crate::coordinates::frames::ECEF;
    use crate::time::{JulianDate, ModifiedJulianDate};
    use qtty::*;

    fn greenwich() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(0.0 * DEG, 51.48 * DEG, 0.0 * M)
    }

    fn one_month() -> Period<MJD> {
        let start = ModifiedJulianDate::from(JulianDate::J2000);
        Period::new(start, start + Days::new(30.0))
    }

    #[test]
    fn phase_topocentric_illuminated_fraction_bounded() {
        let geom = Moon::phase_topocentric(JulianDate::J2000, greenwich());
        assert!(geom.illuminated_fraction >= 0.0 && geom.illuminated_fraction <= 1.0);
    }

    #[test]
    fn illumination_above_returns_periods() {
        // Any fraction above 0 must find at least some time windows in 30 days
        let periods = Moon::illumination_above(one_month(), 0.0, PhaseSearchOpts::default());
        // 0% minimum means always above (all illumination ≥ 0)
        assert!(!periods.is_empty());
    }

    #[test]
    fn illumination_below_returns_periods() {
        // Any fraction below 1.0 must find at least some time windows in 30 days
        let periods = Moon::illumination_below(one_month(), 1.0, PhaseSearchOpts::default());
        assert!(!periods.is_empty());
    }

    #[test]
    fn illumination_range_returns_periods() {
        let periods = Moon::illumination_range(one_month(), 0.0, 1.0, PhaseSearchOpts::default());
        // [0.0, 1.0] covers the entire range — should cover the full window
        assert!(!periods.is_empty());
    }

    #[test]
    fn illumination_above_empty_when_impossible() {
        // k_min above 1.0 — no illumination can exceed 100%
        let periods = Moon::illumination_above(one_month(), 1.01, PhaseSearchOpts::default());
        assert!(periods.is_empty());
    }

    #[test]
    fn illumination_below_empty_when_impossible() {
        // k_max below 0.0 — illumination is never negative
        let periods = Moon::illumination_below(one_month(), -0.01, PhaseSearchOpts::default());
        assert!(periods.is_empty());
    }
}
