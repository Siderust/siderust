// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar apparent-position equations
//!
//! ## Scientific scope
//!
//! The Sun's *apparent* position as seen from Earth differs from its
//! geometric/mean position by three main effects: (1) precession of the
//! equinoxes (continuous drift of the equatorial frame), (2) nutation
//! (periodic perturbation of Earth's spin axis due to the Moon and Sun),
//! and (3) light aberration (apparent displacement caused by Earth's orbital
//! velocity, ~20.5 arc-seconds). This module computes a subset of those
//! corrections that is sufficient for solar position, sunrise/sunset, and
//! shadow modeling: the geometric VSOP87 position is corrected for precession
//! and nutation (via the IAU 2006/2000B NPB matrix), but aberration is
//! intentionally omitted. For topocentric coordinates the full pipeline
//! includes the diurnal parallax correction (~8.9 arc-seconds at the horizon).
//!
//! ## Technical scope
//!
//! Methods on the [`Sun`](crate::bodies::solar_system::Sun) zero-size marker:
//!
//! - [`Sun::get_apparent_geocentric_equ`] — geocentric equatorial apparent
//!   position (`EquatorialTrueOfDate`), nutation applied.
//! - [`Sun::get_apparent_topocentric_equ`] — topocentric apparent position,
//!   adding diurnal parallax shift for a given [`Geodetic<ECEF>`] site.
//! - [`Sun::get_horizontal`] — altitude/azimuth as seen from a site; accepts
//!   any type convertible to `JulianDate` (e.g. `ModifiedJulianDate`).
//! - [`Sun::ecliptic_longitude_geocentric`] — raw geocentric ecliptic
//!   longitude in the J2000 mean ecliptic frame, returned as [`Radians`].
//!
//! All methods are generic over the length unit `U` (pass
//! `AstronomicalUnit`, `Meter`, etc.). The VSOP87 barycentric ephemeris
//! is evaluated via the [`Transform`](crate::coordinates::transform::Transform)
//! blanket impl.
//!
//! ## References
//!
//! - Bretagnon, P., & Francou, G. (1988). Planetary theories in rectangular
//!   and spherical variables: VSOP87 solution. *A&A* **202**, 309–315.
//! - Capitaine, N., & Wallace, P. T. (2006). High precision methods for
//!   locating the celestial intermediate pole and origin. *A&A* **450**, 855.
//! - IERS Conventions (2010), §5.5: Nutation matrix.
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Ch. 27. Willmann-Bell.

use crate::bodies::solar_system::Sun;

use crate::astro::nutation::nutation_iau2000b;
use crate::astro::precession;
use crate::coordinates::{
    cartesian,
    centers::*,
    frames, spherical,
    transform::{Transform, TransformFrame},
};
use crate::event::horizontal;
use crate::qtty::{
    AstronomicalUnit, AstronomicalUnits, LengthUnit, Meter, Quantity, Radian, Radians,
};
use crate::time::JulianDate;

impl Sun {
    /// Returns the **apparent geocentric equatorial coordinates** of the Sun
    /// at a given Julian day.
    ///
    /// This method accounts for:
    /// - **Nutation** in longitude (due to lunar/solar perturbations of Earth's axis).
    ///
    /// ### What is **not** included
    /// - **Aberration** of light: the ~20.5″ displacement caused by Earth's
    ///   orbital velocity is not subtracted.  For most applications the
    ///   geometric direction is sufficient; add aberration separately if needed.
    /// - **Light-time** corrections and relativistic effects.
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
        let helio = cartesian::position::EclipticMeanJ2000::<U, Heliocentric>::CENTER;
        // Center shift first (stays in EclipticMeanJ2000), then frame-rotate to equatorial.
        // Doing both in one `transform` call would apply the ecliptic shift vector into the
        // already-rotated equatorial frame, yielding a near-zero Dec instead of ~-23°.
        let geo_ecl: cartesian::position::EclipticMeanJ2000<U, Geocentric> = helio.transform(jd);
        let geo_cart: cartesian::position::EquatorialMeanJ2000<U, Geocentric> = geo_ecl.to_frame();

        // Apply full IAU 2006/2000B NPB matrix (GCRS → true equator/equinox of date)
        let nut = nutation_iau2000b(jd);
        let npb = precession::precession_nutation_matrix(jd, nut.dpsi, nut.deps);
        let [x_t, y_t, z_t] = npb * [geo_cart.x(), geo_cart.y(), geo_cart.z()];

        let true_cart =
            cartesian::Position::<Geocentric, frames::EquatorialTrueOfDate, U>::new(x_t, y_t, z_t);
        spherical::Position::from_cartesian(&true_cart)
    }

    /// Returns the Sun's apparent topocentric equatorial coordinates as seen
    /// from a given `Geodetic<ECEF>` at the specified Julian Date.
    pub fn get_apparent_topocentric_equ<U: LengthUnit>(
        jd: JulianDate,
        site: Geodetic<frames::ECEF>,
    ) -> spherical::Position<Topocentric, frames::EquatorialTrueOfDate, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    {
        // 1) Compute geocentric cartesian in J2000 (mean) as base.
        // Center shift first (stays in EclipticMeanJ2000), then frame-rotate to equatorial.
        let helio = cartesian::position::EclipticMeanJ2000::<U, Heliocentric>::CENTER;
        let geo_ecl_j2000: cartesian::position::EclipticMeanJ2000<U, Geocentric> =
            helio.transform(jd);
        let geo_cart_j2000: cartesian::position::EquatorialMeanJ2000<U, Geocentric> =
            geo_ecl_j2000.to_frame();

        // 2-5) Shared pipeline: topocentric parallax → precession → nutation → spherical
        horizontal::geocentric_j2000_to_apparent_topocentric(&geo_cart_j2000, site, jd)
    }

    /// Returns the Sun's **horizontal coordinates** (altitude, azimuth) as seen
    /// from a given `Geodetic<ECEF>` at the specified Julian Date or Modified Julian Date.
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
    /// use siderust::coordinates::centers::Geodetic;
    /// use siderust::coordinates::frames::ECEF;
    /// use siderust::time::{JulianDate, ModifiedJulianDate};
    /// use siderust::qtty::*;
    ///
    /// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
    ///
    /// // Using JulianDate
    /// let sun_pos = Sun::get_horizontal::<AstronomicalUnit>(siderust::J2000, site);
    /// println!("Sun altitude: {}", sun_pos.alt().to::<Deg>());
    ///
    /// // Using ModifiedJulianDate
    /// let mjd = siderust::ModifiedJulianDate::new(60000.0);
    /// let sun_pos =
    ///     Sun::get_horizontal::<AstronomicalUnit>(mjd.to::<siderust::time::JD>(), site);
    /// ```
    pub fn get_horizontal<U: LengthUnit>(
        time: impl Into<JulianDate>,
        site: Geodetic<frames::ECEF>,
    ) -> spherical::Position<Topocentric, frames::Horizontal, U>
    where
        Quantity<U>: From<Quantity<Meter>> + From<AstronomicalUnits>,
    {
        let jd = time.into();
        let eq = Self::get_apparent_topocentric_equ::<U>(jd, site);
        horizontal::equatorial_to_horizontal(&eq, site, jd)
    }

    /// Geocentric apparent ecliptic longitude of the Sun in the J2000 mean
    /// ecliptic frame, computed via VSOP87.
    ///
    /// The Sun is placed at the solar-system barycentre / heliocentric origin,
    /// then transformed to a geocentric observer using the VSOP87 cartesian
    /// ephemeris for Earth.  The resulting ecliptic cartesian vector is
    /// converted to spherical and the longitude component is returned.
    ///
    /// ### Returns
    /// Longitude in [`Radians`] (`Quantity<Radian>`).  Convert to degrees with
    /// `.to::<Deg>()` or any other angular unit.
    ///
    /// ### Reference
    /// VSOP87 (Bretagnon & Francou 1988); frame: `EclipticMeanJ2000`.
    pub fn ecliptic_longitude_geocentric(jd: JulianDate) -> Radians {
        let helio =
            cartesian::position::EclipticMeanJ2000::<AstronomicalUnit, Heliocentric>::CENTER;
        let geo_ecl: cartesian::position::EclipticMeanJ2000<AstronomicalUnit, Geocentric> =
            helio.transform(jd);
        spherical::Position::from_cartesian(&geo_ecl)
            .direction()
            .lon()
            .to::<Radian>()
    }
}

#[cfg(test)]
mod tests {
    use crate::bodies::solar_system::Sun;
    use crate::qtty::{AstronomicalUnit, AstronomicalUnits, Degrees, Radians};

    #[test]
    fn apparent_sun_position_j2000() {
        let pos = Sun::get_apparent_geocentric_equ::<AstronomicalUnit>(crate::J2000);

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

    /// At J2000.0 the geocentric ecliptic longitude of the Sun is roughly
    /// 280.46° ≈ 4.895 rad (Sun near winter solstice as seen from Earth).
    #[test]
    fn ecliptic_longitude_geocentric_j2000() {
        let lon: Radians = Sun::ecliptic_longitude_geocentric(crate::J2000);
        // VSOP87 gives ~4.8935 rad (≈ 280.34°) at J2000.0.  Tolerance: 1e-3 rad (~0.06°).
        let expected = 4.8935_f64;
        assert!(
            (lon.value() - expected).abs() < 1e-3,
            "ecliptic longitude mismatch: got {:.6} rad, expected ~{:.4} rad",
            lon.value(),
            expected
        );
    }
}
