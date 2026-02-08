// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::Sun;

use crate::astro::nutation::nutation_rotation;
use crate::astro::{nutation::corrected_ra_with_nutation, precession, JulianDate};
use crate::coordinates::transform::centers::ToTopocentricExt;
use crate::coordinates::{cartesian, centers::*, frames, spherical, transform::Transform};
use qtty::{AstronomicalUnits, LengthUnit, Meter, Quantity};

impl Sun {
    /// Returns the **apparent geocentric equatorial coordinates** of the Sun
    /// at a given Julian day.
    ///
    /// This method accounts for:
    /// - **Nutation** in longitude (due to lunar/solar perturbations of Earth's axis).
    /// - **Aberration** of light caused by Earth's orbital velocity.
    ///
    /// ### Parameters
    /// - `jd`: Julian Day for which to compute the Sun’s apparent position.
    ///
    /// ### Returns
    /// - A `spherical::Position<Geocentric, EquatorialTrueOfDate>` representing the Sun’s
    ///   apparent right ascension and declination, in degrees.
    ///
    /// ### Notes
    /// - This is a simplified model:
    ///   - The heliocentric position of the Sun is treated as the origin.
    ///   - Light-time corrections and relativistic effects are not applied.
    ///   - Nutation and aberration are applied as scalar corrections to the
    ///     azimuthal (longitude) coordinate.
    ///
    /// ### Accuracy
    /// Suitable for applications where approximate solar position is acceptable,
    /// such as sunrise/sunset estimation, shadow modeling, or general astronomy
    /// visualization.
    pub fn get_apparent_geocentric_equ<U: LengthUnit>(
        jd: JulianDate,
    ) -> spherical::position::EquatorialTrueOfDate<U>
    where
        Quantity<U>: From<AstronomicalUnits>,
    {
        let helio = cartesian::position::Ecliptic::<U, Heliocentric>::CENTER;
        let geo_cart: cartesian::position::EquatorialMeanJ2000<U, Geocentric> = helio.transform(jd);
        let geo_sph = spherical::Position::from_cartesian(&geo_cart);
        let mean_of_date = precession::precess_from_j2000(geo_sph, jd);
        let ra = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
        spherical::position::EquatorialTrueOfDate::<U>::new(
            ra,
            mean_of_date.dec(),
            mean_of_date.distance(),
        )
    }

    /// Returns the Sun's apparent topocentric equatorial coordinates as seen
    /// from a given `ObserverSite` at the specified Julian Date.
    pub fn get_apparent_topocentric_equ<U: LengthUnit>(
        jd: JulianDate,
        site: ObserverSite,
    ) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    {
        // 1) Compute geocentric cartesian in J2000 (mean) as base
        let helio = cartesian::position::Ecliptic::<U, Heliocentric>::CENTER;
        let geo_cart_j2000: cartesian::position::EquatorialMeanJ2000<U, Geocentric> =
            helio.transform(jd);

        // 2) Translate geocentric -> topocentric in J2000 frame (applies parallax)
        let topo_cart_j2000 = geo_cart_j2000.to_topocentric(site, jd);

        // 3) Rotate J2000 -> mean-of-date using precession rotation
        let rot_prec = precession::precession_rotation_from_j2000(jd);
        let [x_m, y_m, z_m] = rot_prec.apply_array([
            topo_cart_j2000.x().value(),
            topo_cart_j2000.y().value(),
            topo_cart_j2000.z().value(),
        ]);

        let topo_cart_mod =
            cartesian::Position::<Topocentric, frames::EquatorialMeanOfDate, U>::new_with_params(
                *topo_cart_j2000.center_params(),
                Quantity::<U>::new(x_m),
                Quantity::<U>::new(y_m),
                Quantity::<U>::new(z_m),
            );

        // 4) Apply nutation rotation (mean-of-date -> true-of-date)
        let rot_nut = nutation_rotation(jd);
        let [x_t, y_t, z_t] = rot_nut.apply_array([
            topo_cart_mod.x().value(),
            topo_cart_mod.y().value(),
            topo_cart_mod.z().value(),
        ]);

        let topo_cart_true =
            cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, U>::new_with_params(
                *topo_cart_mod.center_params(),
                Quantity::<U>::new(x_t),
                Quantity::<U>::new(y_t),
                Quantity::<U>::new(z_t),
            );

        // 5) Convert to spherical topocentric equatorial (true of date)
        spherical::Position::from_cartesian(&topo_cart_true)
    }

    /// Returns the Sun's **horizontal coordinates** (altitude, azimuth) as seen
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
    /// use siderust::bodies::solar_system::Sun;
    /// use siderust::coordinates::centers::ObserverSite;
    /// use siderust::astro::{JulianDate, ModifiedJulianDate};
    /// use qtty::*;
    ///
    /// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Using JulianDate
    /// let sun_pos = Sun::get_horizontal::<AstronomicalUnit>(JulianDate::J2000, site);
    /// println!("Sun altitude: {:.2}°", sun_pos.alt().value());
    ///
    /// // Using ModifiedJulianDate
    /// let mjd = ModifiedJulianDate::new(60000.0);
    /// let sun_pos = Sun::get_horizontal::<AstronomicalUnit>(mjd, site);
    /// ```
    pub fn get_horizontal<U: LengthUnit>(
        time: impl Into<JulianDate>,
        site: ObserverSite,
    ) -> spherical::Position<Topocentric, frames::Horizontal, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    {
        use crate::astro::sidereal::{calculate_gst, calculate_lst};
        use qtty::Radian;

        let jd = time.into();
        let sun_eq = Self::get_apparent_topocentric_equ::<U>(jd, site);

        // Get RA and Dec from spherical position (using generic accessors for Topocentric)
        let ra = sun_eq.azimuth();
        let dec = sun_eq.polar();
        let distance = sun_eq.distance();

        // Compute hour angle
        let gst = calculate_gst(jd);
        let lst = calculate_lst(gst, site.lon);
        let ha = (lst - ra).normalize().to::<Radian>();

        // Convert equatorial to horizontal
        let lat = site.lat.to::<Radian>();
        let dec_rad = dec.to::<Radian>();

        // Altitude: sin(alt) = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(HA)
        let sin_alt = dec_rad.sin() * lat.sin() + dec_rad.cos() * lat.cos() * ha.cos();
        let alt = qtty::Degrees::new(sin_alt.asin().to_degrees());

        // Azimuth: tan(az) = -sin(HA) / (cos(lat)*tan(dec) - sin(lat)*cos(HA))
        let az_rad = (-dec_rad.cos() * ha.sin())
            .atan2(dec_rad.sin() * lat.cos() - dec_rad.cos() * ha.cos() * lat.sin());
        let az = qtty::Degrees::new(az_rad.to_degrees()).normalize();

        spherical::Position::<Topocentric, frames::Horizontal, U>::new_with_site(
            site, alt, az, distance,
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::astro::JulianDate;
    use crate::bodies::solar_system::Sun;
    use qtty::AstronomicalUnit;

    #[test]
    fn apparent_sun_position_j2000() {
        let pos = Sun::get_apparent_geocentric_equ::<AstronomicalUnit>(JulianDate::J2000);

        // Expected approximate values around J2000 epoch
        let expected_ra = 281.2; // degrees
        let expected_dec = -23.0; // degrees
        let expected_dist = 1.0; // astronomical units

        eprintln!(
            "Got RA: {}, Dec: {}, Dist: {}",
            pos.ra().value(),
            pos.dec().value(),
            pos.distance.value()
        );

        assert!(
            (pos.ra().value() - expected_ra).abs() < 2.0,
            "RA mismatch: got {}, expected ~{}",
            pos.ra().value(),
            expected_ra
        );
        assert!(
            (pos.dec().value() - expected_dec).abs() < 2.0,
            "Dec mismatch: got {}, expected ~{}",
            pos.dec().value(),
            expected_dec
        );
        assert!((pos.distance.value() - expected_dist).abs() < 0.2);
    }
}
