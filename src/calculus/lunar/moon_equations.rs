// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::Moon;

use crate::astro::nutation::nutation_rotation;
use crate::astro::{precession, JulianDate};
use crate::coordinates::transform::centers::ToTopocentricExt;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{cartesian, centers::*, frames, spherical};
use qtty::{Degree, Kilometer, Kilometers, LengthUnit, Meter, Quantity, Radians};

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
        Quantity<U>: From<Quantity<Meter>> + From<Quantity<Kilometer>>,
    {
        // 1) Get Moon's geocentric ecliptic position from ELP2000
        let moon_geo_ecliptic: cartesian::Position<Geocentric, frames::Ecliptic, Kilometer> =
            Moon::get_geo_position(jd);

        // 2) Transform: Ecliptic → EquatorialMeanJ2000
        let moon_geo_eq_j2000: cartesian::Position<
            Geocentric,
            frames::EquatorialMeanJ2000,
            Kilometer,
        > = TransformFrame::to_frame(&moon_geo_ecliptic);

        // 3) Apply topocentric parallax correction (critical for Moon!)
        let moon_topo_eq_j2000 = moon_geo_eq_j2000.to_topocentric(site, jd);

        // 4) Apply precession: J2000 → mean-of-date
        let rot_prec = precession::precession_rotation_from_j2000(jd);
        let [x_m, y_m, z_m] = rot_prec.apply_array([
            moon_topo_eq_j2000.x().value(),
            moon_topo_eq_j2000.y().value(),
            moon_topo_eq_j2000.z().value(),
        ]);

        let moon_topo_mod = cartesian::Position::<
            Topocentric,
            frames::EquatorialMeanOfDate,
            Kilometer,
        >::new_with_params(
            *moon_topo_eq_j2000.center_params(),
            Kilometers::new(x_m),
            Kilometers::new(y_m),
            Kilometers::new(z_m),
        );

        // 5) Apply nutation rotation (mean-of-date → true-of-date)
        let rot_nut = nutation_rotation(jd);
        let [x_t, y_t, z_t] = rot_nut.apply_array([
            moon_topo_mod.x().value(),
            moon_topo_mod.y().value(),
            moon_topo_mod.z().value(),
        ]);

        let moon_topo_true = cartesian::Position::<
            Topocentric,
            frames::EquatorialTrueOfDate,
            Kilometer,
        >::new_with_params(
            *moon_topo_mod.center_params(),
            Kilometers::new(x_t),
            Kilometers::new(y_t),
            Kilometers::new(z_t),
        );

        // 6) Convert to target unit U and return spherical position
        let x_u: Quantity<U> = moon_topo_true.x().into();
        let y_u: Quantity<U> = moon_topo_true.y().into();
        let z_u: Quantity<U> = moon_topo_true.z().into();

        let moon_topo_true_u =
            cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, U>::new_with_params(
                *moon_topo_true.center_params(),
                x_u,
                y_u,
                z_u,
            );

        spherical::Position::from_cartesian(&moon_topo_true_u)
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
    /// use siderust::astro::{JulianDate, ModifiedJulianDate};
    /// use qtty::*;
    ///
    /// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Using JulianDate
    /// let moon_pos = Moon::get_horizontal::<Kilometer>(JulianDate::J2000, site);
    /// println!("Moon altitude: {:.2}°", moon_pos.alt().value());
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
        Quantity<U>: From<Quantity<Meter>> + From<Quantity<Kilometer>>,
    {
        use crate::astro::sidereal::{calculate_gst, calculate_lst};
        use qtty::Radian;

        let jd = time.into();
        let moon_eq = Self::get_apparent_topocentric_equ::<U>(jd, site);

        // Get RA and Dec from spherical position
        let ra = moon_eq.azimuth();
        let dec = moon_eq.polar();
        let distance = moon_eq.distance();

        // Compute hour angle
        let gst = calculate_gst(jd);
        let lst = calculate_lst(gst, site.lon);
        let ha = (lst - ra).normalize();

        // Convert equatorial to horizontal using qtty units
        let lat_rad = site.lat.to::<Radian>();
        let dec_rad = dec.to::<Radian>();
        let ha_rad = ha.to::<Radian>();

        // Altitude: sin(alt) = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(HA)
        let sin_alt = dec_rad.sin() * lat_rad.sin() + dec_rad.cos() * lat_rad.cos() * ha_rad.cos();
        let alt = Radians::new(sin_alt.asin()).to::<Degree>();

        // Azimuth: computed from hour angle and declination
        let az_rad = (-dec_rad.cos() * ha_rad.sin())
            .atan2(dec_rad.sin() * lat_rad.cos() - dec_rad.cos() * ha_rad.cos() * lat_rad.sin());
        let az = Radians::new(az_rad).normalize().to::<Degree>();

        spherical::Position::<Topocentric, frames::Horizontal, U>::new_with_site(
            site, alt, az, distance,
        )
    }
}
