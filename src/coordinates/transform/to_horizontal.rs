use crate::astro::sidereal::*;
use crate::astro::JulianDate;
use crate::coordinates::{cartesian, spherical};
use qtty::*;

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
pub fn geocentric_to_horizontal(
    target: &spherical::direction::Equatorial,
    observer: &spherical::position::Geographic,
    jd: JulianDate,
) -> spherical::direction::Horizontal {
    let gst = calculate_gst(jd);
    let lst = calculate_lst(gst, observer.lon());

    let ra_deg = target.ra().normalize();
    let dec_rad = target.dec().to::<Radian>();
    let ha_rad = (lst - ra_deg).normalize().to::<Radian>();
    let lat_rad = observer.lat().to::<Radian>();

    let alt_rad =
        (dec_rad.sin() * lat_rad.sin() + dec_rad.cos() * lat_rad.cos() * ha_rad.cos()).asin();

    let az_rad = (-dec_rad.cos() * ha_rad.sin())
        .atan2(dec_rad.sin() * lat_rad.cos() - dec_rad.cos() * ha_rad.cos() * lat_rad.sin());

    spherical::direction::Horizontal::new(
        Degrees::new(alt_rad.to_degrees()),
        Degrees::new(az_rad.to_degrees()),
    )
}

impl spherical::direction::Equatorial {
    pub fn to_horizontal(
        &self,
        observer: &spherical::position::Geographic,
        jd: JulianDate,
    ) -> spherical::direction::Horizontal {
        geocentric_to_horizontal(self, observer, jd)
    }
}

impl cartesian::direction::Equatorial {
    pub fn to_horizontal(
        &self,
        observer: &spherical::position::Geographic,
        jd: JulianDate,
    ) -> cartesian::direction::Horizontal {
        self.to_spherical()
            .to_horizontal(observer, jd)
            .to_cartesian()
    }
}

impl<U: LengthUnit> spherical::position::Equatorial<U> {
    pub fn to_horizontal(
        &self,
        observer: &spherical::position::Geographic,
        jd: JulianDate,
    ) -> spherical::position::Horizontal<U> {
        self.direction()
            .to_horizontal(observer, jd)
            .position(self.distance)
    }
}

impl<U: LengthUnit> cartesian::position::Equatorial<U> {
    pub fn to_horizontal(
        &self,
        observer: &spherical::position::Geographic,
        jd: JulianDate,
    ) -> cartesian::position::Horizontal<U> {
        self.to_spherical()
            .to_horizontal(observer, jd)
            .to_cartesian()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::macros::assert_spherical_eq;

    #[test]
    fn test_sirius_to_horizontal() {
        use crate::bodies::catalog::SIRIUS;
        use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

        let jd: JulianDate = JulianDate::new(2460677.04358);

        let expected_horizontal = spherical::position::Horizontal::<LightYear>::new(
            Degrees::from_dms(-77, 59, 0.0).to::<Deg>(),
            Degrees::from_dms(349, 24, 0.0).to::<Deg>(),
            SIRIUS.target.get_position().distance,
        );

        let horizontal = SIRIUS
            .target
            .get_position()
            .to_horizontal(&ROQUE_DE_LOS_MUCHACHOS, jd);
        assert_spherical_eq!(horizontal, expected_horizontal, 1.5);
    }
}
