use crate::astro::sidereal::*;
use crate::astro::JulianDate;
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::{cartesian, spherical};
use qtty::*;

/// Converts geocentric equatorial coordinates to topocentric horizontal coordinates
/// for a specific observer at a given Julian Date.
///
/// # Deprecated
///
/// This function is deprecated. Use the `TransformToHorizontal` trait instead:
///
/// ```ignore
/// use siderust::coordinates::transform::TransformToHorizontal;
/// 
/// let topo_eq = equatorial.to_topocentric(site, jd);
/// let horizontal = topo_eq.to_horizontal(jd);
/// ```
///
/// # Parameters
/// - `target`: The celestial target's position in geocentric equatorial coordinates (RA/Dec).
/// - `observer`: The geographic coordinates (latitude, longitude) of the observer.
/// - `jd`: The Julian Date of observation.
///
/// # Returns
/// The topocentric horizontal coordinates of the target (Altitude/Azimuth).
/// The returned coordinate embeds the observer's site (converted from the geographic position).
///
/// Azimuth is measured from North (0°) clockwise toward East (90°–360°).
///
/// # See Also
/// - [`calculate_gst`]
/// - [`calculate_lst`]
#[deprecated(
    since = "0.4.0",
    note = "Use TransformToHorizontal trait instead. See module docs for migration guide."
)]
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

    // Create ObserverSite from geographic position
    // Note: Geographic uses altitude from Earth center in km, ObserverSite uses height above ellipsoid in m
    // This is a simplified conversion; for precise work, additional corrections may be needed
    use crate::bodies::EARTH;
    let height_m = (observer.distance.value() - EARTH.radius.value()) * 1000.0;
    let site = ObserverSite::new(observer.lon(), observer.lat(), Quantity::<Meter>::new(height_m));

    spherical::direction::Horizontal::with_site(
        site,
        Degrees::new(alt_rad.to_degrees()),
        Degrees::new(az_rad.to_degrees()),
    )
}

impl spherical::direction::Equatorial {
    /// Converts geocentric equatorial direction to topocentric horizontal.
    ///
    /// # Deprecated
    ///
    /// Use `TransformToHorizontal` trait instead for type-safe transforms.
    #[deprecated(
        since = "0.4.0",
        note = "Use TransformToHorizontal trait. Convert to Topocentric center first."
    )]
    #[allow(deprecated)]
    pub fn to_horizontal(
        &self,
        observer: &spherical::position::Geographic,
        jd: JulianDate,
    ) -> spherical::direction::Horizontal {
        geocentric_to_horizontal(self, observer, jd)
    }
}

impl cartesian::direction::Equatorial {
    /// Converts geocentric equatorial direction to topocentric horizontal.
    ///
    /// # Deprecated
    ///
    /// Use `TransformToHorizontal` trait instead for type-safe transforms.
    #[deprecated(
        since = "0.4.0",
        note = "Use TransformToHorizontal trait. Convert to Topocentric center first."
    )]
    #[allow(deprecated)]
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
    /// Converts geocentric equatorial position to topocentric horizontal.
    ///
    /// # Deprecated
    ///
    /// Use `TransformToHorizontal` trait instead for type-safe transforms.
    #[deprecated(
        since = "0.4.0",
        note = "Use TransformToHorizontal trait. Convert to Topocentric center first."
    )]
    #[allow(deprecated)]
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
    /// Converts geocentric equatorial position to topocentric horizontal.
    ///
    /// # Deprecated
    ///
    /// Use `TransformToHorizontal` trait instead for type-safe transforms.
    #[deprecated(
        since = "0.4.0",
        note = "Use TransformToHorizontal trait. Convert to Topocentric center first."
    )]
    #[allow(deprecated)]
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
    use crate::coordinates::centers::ObserverSite;
    use crate::macros::assert_spherical_eq;

    #[test]
    #[allow(deprecated)]
    fn test_sirius_to_horizontal() {
        use crate::bodies::catalog::SIRIUS;
        use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

        let jd: JulianDate = JulianDate::new(2460677.04358);

        // Create ObserverSite from the observatory's Geographic position
        let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

        let expected_horizontal = spherical::position::Horizontal::<LightYear>::with_site(
            site.clone(),
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
