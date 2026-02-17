// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::Moon;

use crate::calculus::ephemeris::Ephemeris;
use crate::calculus::horizontal;
use crate::coordinates::transform::context::DefaultEphemeris;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, Kilometer, LengthUnit, Meter, Quantity};

impl Moon {
    /// Returns the **apparent topocentric equatorial coordinates** of the Moon
    /// as seen from a given `ObserverSite` at the specified Julian Date.
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
        site: ObserverSite,
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
    /// from a given `ObserverSite` at the specified Julian Date or Modified Julian Date.
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
    /// use siderust::coordinates::centers::ObserverSite;
    /// use siderust::time::{JulianDate, ModifiedJulianDate};
    /// use qtty::*;
    ///
    /// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
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
        site: ObserverSite,
    ) -> spherical::Position<Topocentric, frames::Horizontal, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<Quantity<Kilometer>> + From<AstronomicalUnits>,
    {
        let jd = time.into();
        let eq = Self::get_apparent_topocentric_equ::<U>(jd, site);
        horizontal::equatorial_to_horizontal(&eq, site, jd)
    }
}
