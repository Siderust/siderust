use crate::bodies::solar_system::Sun;

use crate::astro::{nutation::corrected_ra_with_nutation, JulianDate};
use crate::coordinates::{cartesian, centers::*, spherical, transform::Transform};
use crate::units::{AstronomicalUnits, LengthUnit, Quantity};

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
    /// - A `spherical::Position<Geocentric, Equatorial>` representing the Sun’s
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
    ) -> spherical::position::Equatorial<U>
    where
        Quantity<U>: From<AstronomicalUnits>,
    {
        let helio = cartesian::position::Ecliptic::<U, Heliocentric>::CENTER;
        let geo_cart: cartesian::position::Equatorial<U, Geocentric> = helio.transform(jd);
        let geo_sph = geo_cart.to_spherical();
        let ra = corrected_ra_with_nutation(&geo_sph.direction(), jd);
        spherical::position::Equatorial::<U>::new(ra, geo_sph.dec(), geo_sph.distance)
    }
}

#[cfg(test)]
mod tests {
    use crate::astro::JulianDate;
    use crate::bodies::solar_system::Sun;
    use crate::units::AstronomicalUnit;

    #[test]
    fn apparent_sun_position_j2000() {
        let pos = Sun::get_apparent_geocentric_equ::<AstronomicalUnit>(JulianDate::J2000);

        // Expected approximate values around J2000 epoch
        let expected_ra = 281.2; // degrees
        let expected_dec = -23.0; // degrees
        let expected_dist = 1.0; // astronomical units

        assert!(
            (pos.ra().value() - expected_ra).abs() < 2.0,
            "RA mismatch: {} vs {}",
            pos.ra().value(),
            expected_ra
        );
        assert!(
            (pos.dec().value() - expected_dec).abs() < 2.0,
            "Dec mismatch: {} vs {}",
            pos.dec().value(),
            expected_dec
        );
        assert!(
            (pos.distance.value() - expected_dist).abs() < 0.2,
            "Distance mismatch: {} vs {}",
            pos.distance.value(),
            expected_dist
        );
    }
}
