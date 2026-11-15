// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame transformation: EquatorialMeanOfDate → Horizontal
//!
//! ## Scientific scope
//!
//! The **horizontal (alt-az) frame** is an observer-fixed coordinate system
//! referenced to the local horizon. Altitude is the elevation above the
//! horizon and azimuth is the compass bearing (North-clockwise in siderust's
//! native convention). This frame depends on both the observer's geographic
//! location and the sidereal time at observation.
//!
//! ## Technical scope
//!
//! The transformation from the mean equatorial frame of date uses GMST
//! (Greenwich Mean Sidereal Time) derived from the Julian Date, together with
//! the observer's geodetic latitude and longitude embedded in the topocentric
//! center parameters.
//!
//! For true-of-date equatorial coordinates, see
//! [`crate::coordinates::transform::horizontal`].
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. Ch. 13.
//! - Montenbruck, O. & Pfleger, T. (2000). *Astronomy on the Personal
//!   Computer*. §4.2.

use crate::astro::earth_rotation::gmst_from_tt;
use crate::coordinates::centers::Topocentric;
use crate::coordinates::frames::{EquatorialMeanOfDate, Horizontal};
use crate::coordinates::transform::horizontal::{project_cartesian, unproject_cartesian};
use crate::coordinates::transform::Transform;
use crate::coordinates::{cartesian, spherical};
use crate::qtty::LengthUnit;
use crate::time::JulianDate;
use affn::Rotation3;

/// Transform from EquatorialMeanOfDate to Horizontal frame for Topocentric coordinates.
///
/// This transformation requires the Julian Date to compute the local sidereal time.
/// The observer's site information is taken from the coordinate's center params.
impl<U: LengthUnit> Transform<cartesian::Position<Topocentric, Horizontal, U>>
    for cartesian::Position<Topocentric, EquatorialMeanOfDate, U>
{
    fn transform(&self, jd: JulianDate) -> cartesian::Position<Topocentric, Horizontal, U> {
        let site = self.center_params();
        let rotation = Rotation3::rz(-gmst_from_tt(jd));
        let [n, e, u] = project_cartesian(self.as_array().map(|axis| axis.value()), site, rotation);
        cartesian::Position::new_with_params(*site, n, e, u)
    }
}

// =============================================================================
// Horizontal → EquatorialMeanOfDate (for Topocentric center)
// =============================================================================

/// Transform from Horizontal to EquatorialMeanOfDate frame for Topocentric coordinates.
impl<U: LengthUnit> Transform<cartesian::Position<Topocentric, EquatorialMeanOfDate, U>>
    for cartesian::Position<Topocentric, Horizontal, U>
{
    fn transform(
        &self,
        jd: JulianDate,
    ) -> cartesian::Position<Topocentric, EquatorialMeanOfDate, U> {
        let site = self.center_params();
        let rotation = Rotation3::rz(-gmst_from_tt(jd));
        let [x, y, z] =
            unproject_cartesian(self.as_array().map(|axis| axis.value()), site, rotation);
        cartesian::Position::new_with_params(*site, x, y, z)
    }
}

// =============================================================================
// Spherical position implementations
// =============================================================================

impl<U: LengthUnit> Transform<spherical::Position<Topocentric, Horizontal, U>>
    for spherical::Position<Topocentric, EquatorialMeanOfDate, U>
{
    fn transform(&self, jd: JulianDate) -> spherical::Position<Topocentric, Horizontal, U> {
        let cart: cartesian::Position<Topocentric, Horizontal, U> =
            self.to_cartesian().transform(jd);
        spherical::Position::from_cartesian(&cart)
    }
}

impl<U: LengthUnit> Transform<spherical::Position<Topocentric, EquatorialMeanOfDate, U>>
    for spherical::Position<Topocentric, Horizontal, U>
{
    fn transform(
        &self,
        jd: JulianDate,
    ) -> spherical::Position<Topocentric, EquatorialMeanOfDate, U> {
        let cart: cartesian::Position<Topocentric, EquatorialMeanOfDate, U> =
            self.to_cartesian().transform(jd);
        spherical::Position::from_cartesian(&cart)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{centers::Geodetic, frames::ECEF};
    use crate::qtty::{AstronomicalUnit, Degree, Degrees, Radian, AU, DEG, M, RAD};

    type Equatorial = spherical::Position<Topocentric, EquatorialMeanOfDate, AstronomicalUnit>;
    type Horizon = spherical::Position<Topocentric, Horizontal, AstronomicalUnit>;

    fn meridian_zenith(latitude: Degrees) {
        let site = Geodetic::<ECEF>::new(0.0 * DEG, latitude, 0.0 * M);
        let jd = crate::J2000;
        let ra = gmst_from_tt(jd).to::<Degree>();
        let equatorial = Equatorial::new_with_params(site, latitude, ra, 1.0 * AU);
        let horizontal: Horizon = equatorial.transform(jd);
        assert!(
            (horizontal.alt().to::<Radian>() - core::f64::consts::FRAC_PI_2 * RAD).abs()
                < 1e-6 * RAD
        );
    }

    #[test]
    fn test_equatorial_to_horizontal_zenith() {
        meridian_zenith(45.0 * DEG);
    }

    #[test]
    fn test_roundtrip_angles() {
        let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
        let jd = JulianDate::new(2_460_677.0);
        let equatorial = Equatorial::new_with_params(site, -16.72 * DEG, 101.29 * DEG, 1.0 * AU);
        let horizontal: Horizon = equatorial.transform(jd);
        let back: Equatorial = horizontal.transform(jd);
        assert!(
            back.direction()
                .angular_separation(&equatorial.direction())
                .to::<Radian>()
                < 1e-6 * RAD
        );
    }

    // =====================================================================
    // Cartesian position transform roundtrips
    // =====================================================================

    #[test]
    fn test_cartesian_eq_to_horizontal_roundtrip() {
        use crate::coordinates::transform::Transform;
        use crate::qtty::AU;

        let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
        let jd = crate::time::JulianDate::new(2460677.0);

        let eq_pos = cartesian::Position::<
            Topocentric,
            EquatorialMeanOfDate,
            crate::qtty::AstronomicalUnit,
        >::from_array(site, [0.5 * AU, 0.3 * AU, 0.7 * AU]);

        let hz_pos: cartesian::Position<Topocentric, Horizontal, crate::qtty::AstronomicalUnit> =
            eq_pos.transform(jd);
        let eq_back: cartesian::Position<
            Topocentric,
            EquatorialMeanOfDate,
            crate::qtty::AstronomicalUnit,
        > = hz_pos.transform(jd);

        // Distance should be preserved
        assert!(
            (eq_pos.distance() - hz_pos.distance()).abs() < 1e-10 * AU,
            "Distance should be preserved in eq→hz"
        );
        assert!(
            (eq_pos.distance() - eq_back.distance()).abs() < 1e-10 * AU,
            "Distance should be preserved in roundtrip"
        );

        // Position should roundtrip
        assert!((eq_back.x() - eq_pos.x()).abs() < 1e-6 * AU, "x mismatch");
        assert!((eq_back.y() - eq_pos.y()).abs() < 1e-6 * AU, "y mismatch");
        assert!((eq_back.z() - eq_pos.z()).abs() < 1e-6 * AU, "z mismatch");
    }

    // =====================================================================
    // Spherical position transforms
    // =====================================================================

    #[test]
    fn test_spherical_eq_to_horizontal_roundtrip() {
        use crate::coordinates::transform::Transform;
        use crate::qtty::AU;

        let site = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
        let jd = crate::time::JulianDate::new(2460677.0);

        // Build a spherical equatorial position via cartesian → spherical
        let eq_cart = cartesian::Position::<
            Topocentric,
            EquatorialMeanOfDate,
            crate::qtty::AstronomicalUnit,
        >::from_array(site, [0.5 * AU, 0.3 * AU, 0.7 * AU]);
        let eq_sph = spherical::Position::<
            Topocentric,
            EquatorialMeanOfDate,
            crate::qtty::AstronomicalUnit,
        >::from_cartesian(&eq_cart);

        let hz_sph: spherical::Position<Topocentric, Horizontal, crate::qtty::AstronomicalUnit> =
            eq_sph.transform(jd);
        let eq_back: spherical::Position<
            Topocentric,
            EquatorialMeanOfDate,
            crate::qtty::AstronomicalUnit,
        > = hz_sph.transform(jd);

        // Distance should be preserved
        assert!(
            (eq_sph.distance - hz_sph.distance).abs() < 1e-10 * AU,
            "Distance should be preserved"
        );

        // Position should roundtrip
        let eq_cart_back = eq_back.to_cartesian();
        let eq_cart_orig = eq_sph.to_cartesian();
        assert!((eq_cart_back.x() - eq_cart_orig.x()).abs() < 1e-6 * AU, "x");
        assert!((eq_cart_back.y() - eq_cart_orig.y()).abs() < 1e-6 * AU, "y");
        assert!((eq_cart_back.z() - eq_cart_orig.z()).abs() < 1e-6 * AU, "z");
    }

    // =====================================================================
    // Different observer latitudes
    // =====================================================================

    #[test]
    fn test_observer_at_equator() {
        meridian_zenith(0.0 * DEG);
    }
}
