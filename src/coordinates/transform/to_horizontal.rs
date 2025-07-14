use crate::astro::sidereal::*;
use crate::coordinates::{
    spherical, cartesian,
    centers::*,  frames::*
};
use crate::units::{Quantity, LengthUnit, Degrees};
use crate::astro::JulianDate;


/// Converts geocentric equatorial coordinates to topocentric horizontal coordinates
/// for a specific observer at a given Julian Date.
///
/// # Parameters
/// - `target`: The celestial target's position in geocentric equatorial coordinates (RA/Dec).
/// - `observer`: The geographic coordinates (latitude, longitude) of the observer.
/// - `jd`: The Julian Date of observation.
///
/// # Returns
/// The topocentric horizontal coordinates of the target (Altitude/Azimuth).
///
/// Azimuth is measured from North (0°) clockwise toward East (90°–360°).
///
/// # See Also
/// - [`calculate_gst`]
/// - [`calculate_lst`]
pub fn geocentric_to_horizontal<U: LengthUnit>(
    target:   &spherical::position::Equatorial<U>,
    observer: &spherical::position::Geographic,
    jd:       JulianDate
) -> spherical::position::Horizontal<U> {

    // 2) Tiempo sidéreo con ese JD, no con target.t
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, observer.lon());

    // 3) El resto no cambia
    let ra_deg   = target.ra().normalize();
    let dec_rad  = target.dec().to_radians();
    let ha_rad   = (lst - ra_deg).normalize().to_radians();
    let lat_rad  = observer.lat().to_radians();

    let alt_rad = (dec_rad.sin() * lat_rad.sin()
                 + dec_rad.cos() * lat_rad.cos() * ha_rad.cos()).asin();

    // azimut de 0-360°
    let az_rad  = (-dec_rad.cos() * ha_rad.sin()).atan2(
                    dec_rad.sin() * lat_rad.cos()
                  - dec_rad.cos() * ha_rad.cos() * lat_rad.sin());

    spherical::position::Horizontal::new::<Quantity<U>>(
        Degrees::new(alt_rad.to_degrees()),
        Degrees::new(az_rad.to_degrees()),
        target.distance.unwrap(),
    )
}

impl<U: LengthUnit> cartesian::Position<Geocentric, Equatorial, U>
where
    cartesian::Position<Topocentric, Horizontal, U>: for<'a> From<&'a spherical::Position<Topocentric, Horizontal, U>>,
{
    pub fn to_horizontal(&self, observer: &spherical::position::Geographic, jd: JulianDate) -> cartesian::Position<Topocentric, Horizontal, U> {
        let spherical: spherical::Position<Geocentric, Equatorial, U>   = self.into();
        let horizontal = geocentric_to_horizontal(&spherical, observer, jd);
        (&horizontal).into()
    }
}


impl<U: LengthUnit> spherical::Position<Geocentric, Equatorial, U> {
    pub fn to_horizontal(&self, observer: &spherical::position::Geographic, jd: JulianDate) -> spherical::Position<Topocentric, Horizontal, U> {
        geocentric_to_horizontal(self, observer, jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{DMS, LightYear};
    use crate::macros::assert_spherical_eq;

    #[test]
    fn test_sirius_to_horizontal() {
        use crate::bodies::catalog::SIRIUS;
        use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

        let jd: JulianDate = JulianDate::new(2460677.04358);

        let expected_horizontal = spherical::position::Horizontal::<LightYear>::new(
            DMS::new(DMS::NEGATIVE, 77, 59, 0.0).to_degrees(),
            DMS::new(DMS::POSITIVE, 349, 24, 0.0).to_degrees(),
            SIRIUS.target.get_position().distance.unwrap(),
        );

        let horizontal = geocentric_to_horizontal(&SIRIUS.target.get_position(), &ROQUE_DE_LOS_MUCHACHOS, jd);
        assert_spherical_eq!(horizontal, expected_horizontal, 1.5);
    }
}
