use crate::bodies::solar_system::Sun;

use crate::coordinates::{
    cartesian, spherical,
    centers::Geocentric,
    frames::Ecliptic
};
use crate::units::{Quantity, AstronomicalUnits, LengthUnit, Degrees};
use crate::astro::JulianDate;

use crate::astro::nutation::get_nutation;
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
        Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    {
        let helio = cartesian::position::Ecliptic::<U>::CENTER;
        let geo_cart: cartesian::Position<Geocentric, Ecliptic, U> = (&helio).into();
        let mut geo = geo_cart.to_spherical();

        // Apply nutation in ecliptic longitude
        let nutation = get_nutation(jd);
        geo.azimuth += nutation.longitude;

        // Apply aberration correction (simplified constant formula)
        let aberration = (20.4898 / (360.0 * 60.0 * 60.0)) / geo.distance.value();
        geo.azimuth -= Degrees::new(aberration);

        (&geo).into()
    }
}
