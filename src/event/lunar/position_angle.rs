// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Bright-Limb Position Angle
//!
//! ## Scientific scope
//!
//! The bright-limb position angle points from the Moon's disc centre towards
//! the Sun in the sky plane. It is measured from celestial north through east,
//! in the true equator and equinox of the observation date, in `[0, 2π)`.
//! The origin is either Earth's centre or the supplied terrestrial observer.
//! Epochs are TT; both bodies are evaluated at the same epoch with backend `E`.
//!
//! These are geometric directions rotated by the existing IAU 2006 precession
//! and IAU 2006A nutation/frame-bias providers. The topocentric path subtracts
//! the observer position from **both** bodies, using the existing Earth rotation
//! and default IERS EOP pipeline (including its fallback when data is absent).
//! No light-time, annual/diurnal aberration, gravitational deflection, or
//! atmospheric refraction is added; these are not fully apparent directions.
//! Near conjunction/opposition the orientation is sensitive to position errors.
//!
//! ## Technical scope
//!
//! The two backend-generic functions return `Option<qtty::Radians>` and delegate
//! all bearing geometry and singularity detection to
//! [`affn::spherical::Direction::position_angle_to`]. Coincident or antipodal
//! directions, an undefined lunar north reference at a celestial pole, or a
//! vanishing numerical tangent return `None` according to affn's policy.
//! [`Moon`](crate::bodies::solar_system::Moon) provides default-backend wrappers.
//!
//! The greatest defect of illumination is opposite the bright limb. Derive it
//! from an already computed angle without evaluating the ephemeris again:
//!
//! ```rust
//! use siderust::{moon_bright_limb_position_angle_geocentric, ephemeris::Vsop87Ephemeris};
//! use siderust::qtty::Radians;
//!
//! let bright = moon_bright_limb_position_angle_geocentric::<Vsop87Ephemeris>(siderust::J2000);
//! let defect = bright.map(|angle| (angle + Radians::HALF_TURN).normalize());
//! assert!(defect.is_some());
//! ```
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Chapter 41.
//!   Willmann–Bell.
//! - IERS (2010). *IERS Conventions*, Technical Note 36, Chapter 5.
//! - Siderust (2026). affn spherical bearing primitive:
//!   <https://github.com/Siderust/affn/issues/34>.

use crate::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use crate::coordinates::frames;
use crate::coordinates::transform::centers::position::to_topocentric::to_topocentric_with;
use crate::coordinates::transform::ext::PositionAstroExt;
use crate::coordinates::transform::{AstroContext, TransformFrame};
use crate::coordinates::{cartesian, spherical};
use crate::ephemeris::Ephemeris;
use crate::qtty::{AstronomicalUnit, Kilometer, Radians};
use crate::time::JulianDate;

/// Geocentric bright-limb position angle with compile-time ephemeris backend `E`.
///
/// Both geometric directions use the true equator/equinox of date and Earth's
/// centre. See the [module documentation](self) for corrections and limitations.
///
/// # Arguments
///
/// * `jd`: Observation epoch as a TT Julian Date.
///
/// # Returns
///
/// Angle north through east in `[0, 2π)`, or `None` for undefined geometry
/// as determined by affn (including coincident/antipodal directions and poles).
///
/// # Examples
///
/// ```rust
/// use siderust::{moon_bright_limb_position_angle_geocentric, ephemeris::Vsop87Ephemeris};
/// let angle = moon_bright_limb_position_angle_geocentric::<Vsop87Ephemeris>(siderust::J2000);
/// assert!(angle.unwrap().is_finite());
/// ```
pub fn moon_bright_limb_position_angle_geocentric<E: Ephemeris>(jd: JulianDate) -> Option<Radians> {
    let ctx = AstroContext::<E>::default();
    let (moon, sun) = geocentric_equatorial::<E>(jd);
    let moon_true = moon.to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, &ctx);
    let sun_true = sun.to_frame_with::<frames::EquatorialTrueOfDate, _>(&jd, &ctx);
    spherical::Position::from_cartesian(&moon_true)
        .direction()
        .position_angle_to(&spherical::Position::from_cartesian(&sun_true).direction())
}

/// Topocentric bright-limb position angle with compile-time ephemeris backend `E`.
///
/// Both directions share the observer origin and true equator/equinox of date.
/// Parallax is applied to both the Moon and Sun. See the
/// [module documentation](self) for the correction model and limitations.
///
/// # Arguments
///
/// * `jd`: Observation epoch as a TT Julian Date.
/// * `site`: Observer's geodetic longitude, latitude, and height in ECEF.
///
/// # Returns
///
/// Angle north through east in `[0, 2π)`, or `None` when affn determines
/// that the orientation is undefined. This is relative to celestial north,
/// rather than the local zenith; no parallactic-angle rotation is applied.
///
/// # Examples
///
/// ```rust
/// use siderust::{moon_bright_limb_position_angle_topocentric, ephemeris::Vsop87Ephemeris};
/// use siderust::coordinates::{centers::Geodetic, frames::ECEF};
/// use siderust::qtty::{Degrees, Meters};
/// let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
/// let angle = moon_bright_limb_position_angle_topocentric::<Vsop87Ephemeris>(siderust::J2000, site);
/// assert!(angle.unwrap().is_finite());
/// ```
pub fn moon_bright_limb_position_angle_topocentric<E: Ephemeris>(
    jd: JulianDate,
    site: Geodetic<frames::ECEF>,
) -> Option<Radians> {
    let ctx = AstroContext::<E>::default();
    let (moon, sun) = geocentric_equatorial::<E>(jd);
    let moon_topo = to_topocentric_with(&moon, site, jd, &ctx);
    let sun_topo = to_topocentric_with(&sun, site, jd, &ctx);
    let rotation = crate::coordinates::transform::providers::frame_rotation_with::<
        frames::EquatorialMeanJ2000,
        frames::EquatorialTrueOfDate,
        _,
    >(jd, &ctx);
    let moon_true =
        cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, Kilometer>::from_array(
            site,
            rotation * *moon_topo.as_array(),
        );
    let sun_true = cartesian::Position::<Topocentric, frames::EquatorialTrueOfDate, AstronomicalUnit>::from_array(
        site, rotation * *sun_topo.as_array(),
    );
    spherical::Position::from_cartesian(&moon_true)
        .direction()
        .position_angle_to(&spherical::Position::from_cartesian(&sun_true).direction())
}

// Compute each ephemeris state once; negate Earth–Sun to obtain Sun–Earth.
fn geocentric_equatorial<E: Ephemeris>(
    jd: JulianDate,
) -> (
    cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, Kilometer>,
    cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, AstronomicalUnit>,
) {
    let moon = E::moon_geocentric(jd);
    let earth = E::earth_heliocentric(jd);
    let sun = cartesian::Position::<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit>::new(
        -earth.x(),
        -earth.y(),
        -earth.z(),
    );
    (
        TransformFrame::to_frame(&moon),
        TransformFrame::to_frame(&sun),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Moon;
    use crate::coordinates::centers::{Barycentric, Heliocentric};
    use crate::coordinates::transform::context::DefaultEphemeris;
    use crate::ephemeris::AuPerDay;
    use crate::qtty::{Degree, Degrees, Meters};
    use core::sync::atomic::{AtomicUsize, Ordering};

    static MOON_CALLS: AtomicUsize = AtomicUsize::new(0);
    static EARTH_CALLS: AtomicUsize = AtomicUsize::new(0);

    // Known true-of-date configurations, converted to the ephemeris contract
    // using existing inverse rotations. Other ephemeris operations must not run.
    struct Cardinal<const CASE: u8>;

    fn source_position<U: crate::qtty::LengthUnit>(
        jd: JulianDate,
        ra: f64,
        dec: f64,
    ) -> cartesian::Position<Geocentric, frames::EclipticMeanJ2000, U> {
        let true_pos = spherical::Position::<Geocentric, frames::EquatorialTrueOfDate, U>::new(
            Degrees::new(ra),
            Degrees::new(dec),
            crate::qtty::Quantity::<U>::new(1.0),
        )
        .to_cartesian();
        let mean = true_pos.to_frame_with::<frames::EquatorialMeanJ2000, _>(
            &jd,
            &AstroContext::<DefaultEphemeris>::default(),
        );
        TransformFrame::to_frame(&mean)
    }

    impl<const CASE: u8> Ephemeris for Cardinal<CASE> {
        fn moon_geocentric(
            jd: JulianDate,
        ) -> cartesian::Position<Geocentric, frames::EclipticMeanJ2000, Kilometer> {
            MOON_CALLS.fetch_add(1, Ordering::SeqCst);
            source_position(jd, 0.0, 0.0)
        }

        fn earth_heliocentric(
            jd: JulianDate,
        ) -> cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>
        {
            EARTH_CALLS.fetch_add(1, Ordering::SeqCst);
            let (ra, dec) = match CASE {
                0 => (0.0, 30.0),
                1 => (30.0, 0.0),
                2 => (330.0, 0.0),
                _ => (0.0, 0.0),
            };
            let sun = source_position::<AstronomicalUnit>(jd, ra, dec);
            cartesian::Position::new(-sun.x(), -sun.y(), -sun.z())
        }

        fn sun_barycentric(
            _: JulianDate,
        ) -> cartesian::Position<Barycentric, frames::EclipticMeanJ2000, AstronomicalUnit> {
            panic!("unexpected barycentric Sun evaluation")
        }
        fn earth_barycentric(
            _: JulianDate,
        ) -> cartesian::Position<Barycentric, frames::EclipticMeanJ2000, AstronomicalUnit> {
            panic!("unexpected barycentric Earth evaluation")
        }
        fn earth_barycentric_velocity(
            _: JulianDate,
        ) -> cartesian::Velocity<frames::EclipticMeanJ2000, AuPerDay> {
            panic!("unexpected aberration evaluation")
        }
    }

    fn normalized(angle: Option<Radians>) -> Radians {
        let angle = angle.expect("defined bearing");
        assert!(angle.is_finite());
        assert!(angle >= Radians::new(0.0) && angle < Radians::FULL_TURN);
        angle
    }

    #[test]
    fn selected_backend_cardinal_bearings_and_single_evaluations() {
        let jd = crate::J2000;
        // At latitude zero and height minus the WGS84 equatorial radius the
        // observer vector is zero, providing an independent topocentric oracle.
        let centre = Geodetic::<frames::ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
            Meters::new(-6_378_137.0),
        );
        MOON_CALLS.store(0, Ordering::SeqCst);
        EARTH_CALLS.store(0, Ordering::SeqCst);
        let north = moon_bright_limb_position_angle_geocentric::<Cardinal<0>>(jd).unwrap();
        assert!(north.wrap_signed().abs() < Radians::new(1e-12));
        let north_topo =
            moon_bright_limb_position_angle_topocentric::<Cardinal<0>>(jd, centre).unwrap();
        assert!(north_topo.wrap_signed().abs() < Radians::new(1e-12));
        for (geo, topo, expected) in [
            (
                moon_bright_limb_position_angle_geocentric::<Cardinal<1>>(jd),
                moon_bright_limb_position_angle_topocentric::<Cardinal<1>>(jd, centre),
                90.0,
            ),
            (
                moon_bright_limb_position_angle_geocentric::<Cardinal<2>>(jd),
                moon_bright_limb_position_angle_topocentric::<Cardinal<2>>(jd, centre),
                270.0,
            ),
        ] {
            for angle in [geo, topo] {
                assert!(
                    (normalized(angle).to::<Degree>() - Degrees::new(expected)).abs()
                        < Degrees::new(1e-10)
                );
            }
        }
        assert_eq!(
            moon_bright_limb_position_angle_geocentric::<Cardinal<3>>(jd),
            None
        );
        assert_eq!(
            moon_bright_limb_position_angle_topocentric::<Cardinal<3>>(jd, centre),
            None,
        );
        let surface = Geodetic::<frames::ECEF>::new(
            Degrees::new(12.0),
            Degrees::new(48.0),
            Meters::new(500.0),
        );
        let custom_topo = normalized(moon_bright_limb_position_angle_topocentric::<Cardinal<1>>(
            jd, surface,
        ));
        let default_topo = Moon::bright_limb_position_angle_topocentric(jd, surface).unwrap();
        assert_ne!(custom_topo, default_topo);
        assert_eq!(MOON_CALLS.load(Ordering::SeqCst), 9);
        assert_eq!(EARTH_CALLS.load(Ordering::SeqCst), 9);
    }

    #[test]
    fn meeus_independent_reference() {
        // Meeus bright-limb example, independently reproduced by PyMeeus:
        // https://github.com/architest/pymeeus/blob/master/pymeeus/Moon.py
        // 1992-04-12 00:00 TT: chi = 285.0 degrees (rounded to 0.1 degree).
        // Allow 0.1 degree for rounding and the differing ephemeris/corrections.
        let jd = JulianDate::new(2_448_724.5);
        let angle = normalized(moon_bright_limb_position_angle_geocentric::<DefaultEphemeris>(jd));
        assert!((angle.to::<Degree>() - Degrees::new(285.0)).abs() < Degrees::new(0.1));
        let defect = (angle + Radians::HALF_TURN).normalize();
        assert!((defect.to::<Degree>() - Degrees::new(105.0)).abs() < Degrees::new(0.1));
    }

    #[test]
    fn undefined_and_near_singular_geometry_delegated_to_affn() {
        let origin = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
        );
        for (ra, dec) in [(0.0, 0.0), (180.0, 0.0)] {
            let other = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
                Degrees::new(ra),
                Degrees::new(dec),
            );
            assert_eq!(origin.position_angle_to(&other), None);
        }
        let pole = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
            Degrees::new(0.0),
            Degrees::new(90.0),
        );
        assert_eq!(pole.position_angle_to(&origin), None);
        for ra in [1e-9, 180.0 - 1e-9, 180.0 + 1e-9, 360.0 - 1e-9] {
            let other = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
                Degrees::new(ra),
                Degrees::new(0.0),
            );
            let angle = normalized(origin.position_angle_to(&other));
            let expected = if ra < 180.0 { 90.0 } else { 270.0 };
            assert!((angle.to::<Degree>() - Degrees::new(expected)).abs() < Degrees::new(1e-10));
        }
    }

    #[test]
    fn parallax_epochs_and_convenience_wrappers() {
        let sites = [(0.0, 51.48), (-70.4, -24.6), (151.2, -33.9)].map(|(lon, lat)| {
            Geodetic::<frames::ECEF>::new(Degrees::new(lon), Degrees::new(lat), Meters::new(0.0))
        });
        // Include epochs near the 2024 April new Moon and March full Moon.
        for jd in [
            crate::J2000,
            JulianDate::new(2_460_409.25),
            JulianDate::new(2_460_394.8),
        ] {
            let geo =
                normalized(moon_bright_limb_position_angle_geocentric::<DefaultEphemeris>(jd));
            assert_eq!(Some(geo), Moon::bright_limb_position_angle_geocentric(jd));
            let mut topocentric = [Radians::new(0.0); 3];
            for (i, site) in sites.iter().enumerate() {
                let topo = normalized(moon_bright_limb_position_angle_topocentric::<
                    DefaultEphemeris,
                >(jd, *site));
                assert_eq!(
                    Some(topo),
                    Moon::bright_limb_position_angle_topocentric(jd, *site)
                );
                assert!((geo - topo).abs() > Radians::new(1e-6));
                topocentric[i] = topo;
            }
            assert!((topocentric[0] - topocentric[1]).abs() > Radians::new(1e-6));
            assert!((topocentric[1] - topocentric[2]).abs() > Radians::new(1e-6));
        }
    }
}
