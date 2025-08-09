use crate::bodies::solar_system::Sun;

use crate::coordinates::{
    cartesian, spherical,
    centers::*,
    transform::Transform,
};
use crate::units::{Quantity, AstronomicalUnits, LengthUnit};
use crate::astro::{JulianDate, nutation::corrected_ra_with_nutation};

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
    pub fn get_apparent_geocentric_equ<U: LengthUnit>(jd: JulianDate) -> spherical::position::Equatorial<U>
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
