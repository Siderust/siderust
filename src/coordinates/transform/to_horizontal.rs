use crate::astro::sidereal::*;
use crate::coordinates::{
    SphericalCoord, CartesianCoord,
    GeographicCoord,
    centers::*,  frames::*
};
use crate::units::{Degrees, JulianDay};


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
    target:   &SphericalCoord<Geocentric, Equatorial>,
    observer: &GeographicCoord,
    jd:       JulianDay
) -> SphericalCoord<Topocentric, Horizontal> {

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

    SphericalCoord::<Topocentric, Horizontal>::new(
        Degrees::new(alt_rad.to_degrees()),
        Degrees::new(az_rad.to_degrees()),
        target.radial_distance,
    )
}


impl CartesianCoord<Geocentric, Equatorial> {
    pub fn to_horizontal(&self, observer: &GeographicCoord, jd: JulianDay) -> CartesianCoord<Topocentric, Horizontal> {
        let spherical: SphericalCoord<Geocentric, Equatorial>   = self.into();
        let horizontal: SphericalCoord<Topocentric, Horizontal> = spherical.to_horizontal(observer, jd);
        (&horizontal).into()
    }
}


impl SphericalCoord<Geocentric, Equatorial> {
    pub fn to_horizontal(&self, observer: &GeographicCoord, jd: JulianDay) -> SphericalCoord<Topocentric, Horizontal> {
        geocentric_to_horizontal(self, observer, jd)
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::*;
    use crate::units::DMS;

    /// Helper function to compare floating-point values within a small tolerance
    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn coords_approx_eq(a: &SphericalCoord<impl ReferenceCenter, impl ReferenceFrame>,
                        b: &SphericalCoord<impl ReferenceCenter, impl ReferenceFrame>,
                        epsilon: f64) {
        assert!(approx_eq(a.polar.as_f64(), b.polar.as_f64(), epsilon), "Current Polar {:?}; Expected {:?}", a.polar, b.polar);
        assert!(approx_eq(a.azimuth.as_f64(), b.azimuth.as_f64(), epsilon), "Current Azimuth {:?}; Expected {:?}", a.azimuth, b.azimuth);
    }

    #[test]
    fn test_sirius_to_horizontal() {
        use crate::bodies::catalog::SIRIUS;
        use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

        let jd: JulianDay = JulianDay::new(2460677.04358);

        let expected_horizontal = SphericalCoord::<Topocentric, Horizontal>::new(
            DMS::new(DMS::NEGATIVE, 77, 59, 0.0).to_degrees(),
            DMS::new(DMS::POSITIVE, 349, 24, 0.0).to_degrees(),
            SIRIUS.target.get_position().radial_distance,
        );

        let horizontal = geocentric_to_horizontal(&SIRIUS.target.get_position(), &ROQUE_DE_LOS_MUCHACHOS, jd);
        coords_approx_eq(&horizontal, &expected_horizontal, 1.5);
    }
}
