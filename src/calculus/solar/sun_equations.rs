// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::Sun;

use crate::astro::{nutation::corrected_ra_with_nutation, precession};
use crate::calculus::horizontal;
use crate::coordinates::{cartesian, centers::*, frames, spherical, transform::Transform};
use crate::time::JulianDate;
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
            mean_of_date.distance,
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

        // 2-5) Shared pipeline: topocentric parallax → precession → nutation → spherical
        horizontal::geocentric_j2000_to_apparent_topocentric(&geo_cart_j2000, site, jd)
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
    /// use siderust::time::{JulianDate, ModifiedJulianDate};
    /// use qtty::*;
    ///
    /// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Using JulianDate
    /// let sun_pos = Sun::get_horizontal::<AstronomicalUnit>(JulianDate::J2000, site);
    /// println!("Sun altitude: {}", sun_pos.alt().to::<Deg>());
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
        let jd = time.into();
        let eq = Self::get_apparent_topocentric_equ::<U>(jd, site);
        horizontal::equatorial_to_horizontal(&eq, site, jd)
    }
}

#[cfg(test)]
mod tests {
    use crate::bodies::solar_system::Sun;
    use crate::time::JulianDate;
    use qtty::{AstronomicalUnit, AstronomicalUnits, Degrees};

    #[test]
    fn apparent_sun_position_j2000() {
        let pos = Sun::get_apparent_geocentric_equ::<AstronomicalUnit>(JulianDate::J2000);

        // Expected approximate values around J2000 epoch
        let expected_ra = 281.2; // degrees
        let expected_dec = -23.0; // degrees
        let expected_dist = 1.0; // astronomical units

        eprintln!(
            "Got RA: {}, Dec: {}, Dist: {}",
            pos.ra(),
            pos.dec(),
            pos.distance
        );

        assert!(
            (pos.ra() - Degrees::new(expected_ra)).abs() < Degrees::new(2.0),
            "RA mismatch: got {}, expected ~{}",
            pos.ra(),
            expected_ra
        );
        assert!(
            (pos.dec() - Degrees::new(expected_dec)).abs() < Degrees::new(2.0),
            "Dec mismatch: got {}, expected ~{}",
            pos.dec(),
            expected_dec
        );
        assert!(
            (pos.distance - AstronomicalUnits::new(expected_dist)).abs()
                < AstronomicalUnits::new(0.2)
        );
    }
}
