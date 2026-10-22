// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Mars and Saturn apparent physical orientation
//!
//! ## Scientific scope
//!
//! This module evaluates geocentric physical-orientation angles for Mars and
//! Saturn at a TT epoch:
//!
//! - **Mars:** position angle of the north rotational pole, and bright-limb
//!   position angle (Sunward bearing on the sky), both measured from celestial
//!   north through east in the true equator and equinox of date, normalised to
//!   `[0, 2π)` when defined.
//! - **Saturn:** position angle of the projected north pole / ring axis (same
//!   convention), plus signed ring opening angles as seen from Earth (`B`) and
//!   from the Sun (`B′`). Positive opening means the north face of the rings is
//!   presented to that observer (saturnicentric latitude of the observer above
//!   the ring plane).
//!
//! Pole directions come from the body's IAU WGCCRE linear model
//! ([`HasIauRotation`](crate::astro::HasIauRotation)) evaluated at the
//! **observation** epoch. Planet heliocentric states are evaluated at a
//! light-time retarded epoch (two fixed-point iterations). Earth is taken at
//! the observation epoch. Geocentric and heliocentric vectors are formed in
//! mean ecliptic J2000, then rotated to the true equator/equinox of date with
//! the existing IAU 2006 precession and IAU 2006A nutation/frame-bias
//! providers.
//!
//! These are **light-time geometric** orientations in the true-of-date frame.
//! Annual/diurnal aberration, gravitational deflection, atmospheric
//! refraction, and topocentric parallax are **not** applied. Do not label the
//! results as fully apparent places without those corrections. Near
//! conjunction or when the body lies near a celestial pole the sky-plane
//! bearings may be undefined (`None`).
//!
//! TT is used as the working time scale; IAU rotational rates are defined on
//! TDB. The TT−TDB difference is neglected here (well below the angular
//! accuracy of the linear pole model over light-time intervals).
//!
//! ## Technical scope
//!
//! Backend-generic entry points require [`PlanetaryEphemeris`] so Mars/Saturn
//! and Earth come from the same backend. Sky-plane bearings delegate to
//! [`affn::spherical::Direction::position_angle_to`]. Opening angles are
//! `asin` of the north-pole component of the unit vector from the planet toward
//! the observer (Earth or Sun). Singular or non-finite geometry yields
//! `None` for position angles; opening angles remain defined whenever the
//! direction vectors normalise.
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Chapters 42 and 45.
//!   Willmann–Bell.
//! - Archinal, B. A. et al., *Report of the IAU Working Group on Cartographic
//!   Coordinates and Rotational Elements*.
//! - IERS (2010). *IERS Conventions*, Technical Note 36, Chapter 5.
//! - Siderust affn spherical bearing:
//!   <https://github.com/Siderust/affn/issues/34>.

use crate::astro::{HasIauRotation, IauRotationParams};
use crate::bodies::solar_system::{Mars, Saturn};
use crate::coordinates::centers::{Geocentric, Heliocentric};
use crate::coordinates::frames;
use crate::coordinates::transform::ext::PositionAstroExt;
use crate::coordinates::transform::{AstroContext, TransformFrame};
use crate::coordinates::{cartesian, spherical};
use crate::ephemeris::PlanetaryEphemeris;
use crate::qtty::{velocity::AU_PER_DAY_C, AstronomicalUnit, Days, Radians};
use crate::time::JulianDate;

#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};

/// Fixed light-time iterations (Meeus planetary physical ephemerides use two).
const LIGHT_TIME_ITERS: u8 = 2;

/// Geocentric Mars physical orientation at a TT epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MarsOrientation {
    /// Position angle of the north rotational pole, north through east in
    /// `[0, 2π)`, or `None` when affn reports undefined sky-plane geometry.
    pub north_pole_position_angle: Option<Radians>,
    /// Bright-limb position angle (bearing toward the Sun), same convention, or
    /// `None` when undefined. Greatest defect of illumination is
    /// `(bright + Radians::HALF_TURN).normalize()` when `bright` is defined.
    pub bright_limb_position_angle: Option<Radians>,
}

/// Geocentric Saturn ring / pole orientation at a TT epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SaturnRingOrientation {
    /// Position angle of the projected north pole (ring axis), north through
    /// east in `[0, 2π)`, or `None` when undefined.
    pub north_pole_position_angle: Option<Radians>,
    /// Signed ring opening angle as seen from Earth (saturnicentric latitude of
    /// the Earth). Positive ⇒ north face of the rings toward Earth.
    pub ring_opening_earth: Radians,
    /// Signed ring opening angle as seen from the Sun (saturnicentric latitude
    /// of the Sun). Positive ⇒ north face illuminated.
    pub ring_opening_sun: Radians,
}

/// Geocentric Mars orientation with compile-time planetary ephemeris `E`.
///
/// # Arguments
///
/// * `jd` — observation epoch on the TT axis as a Julian Date.
///
/// # Returns
///
/// [`MarsOrientation`] with optional sky-plane bearings. See the
/// [module documentation](self) for light-time, frames, and omitted corrections.
///
/// # Examples
///
/// ```rust
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::mars_orientation;
///
/// let ori = mars_orientation::<Vsop87Ephemeris>(siderust::J2000);
/// assert!(ori.north_pole_position_angle.unwrap().is_finite());
/// assert!(ori.bright_limb_position_angle.unwrap().is_finite());
/// ```
pub fn mars_orientation<E: PlanetaryEphemeris>(jd: JulianDate) -> MarsOrientation {
    let ctx = AstroContext::<E>::default();
    let (planet_true, sun_true, pole_true, _) =
        planet_sun_pole_true_of_date::<E, Mars, _>(&jd, &ctx, E::mars_heliocentric);
    MarsOrientation {
        north_pole_position_angle: bearing_to_pole(&planet_true, &pole_true),
        bright_limb_position_angle: bearing_to_position(&planet_true, &sun_true),
    }
}

/// Geocentric Saturn ring orientation with compile-time planetary ephemeris `E`.
///
/// # Arguments
///
/// * `jd` — observation epoch on the TT axis as a Julian Date.
///
/// # Returns
///
/// [`SaturnRingOrientation`]. Opening angles are always finite for non-zero
/// planet vectors; the pole position angle may be `None`. See the
/// [module documentation](self) for conventions and corrections.
///
/// # Examples
///
/// ```rust
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::saturn_ring_orientation;
/// use siderust::qtty::Degrees;
///
/// let ori = saturn_ring_orientation::<Vsop87Ephemeris>(siderust::J2000);
/// assert!(ori.north_pole_position_angle.unwrap().is_finite());
/// assert!(ori.ring_opening_earth.to::<siderust::qtty::Degree>().abs() < Degrees::new(90.0));
/// ```
pub fn saturn_ring_orientation<E: PlanetaryEphemeris>(jd: JulianDate) -> SaturnRingOrientation {
    let ctx = AstroContext::<E>::default();
    let (planet_true, _sun_true, pole_true, planet_helio_true) =
        planet_sun_pole_true_of_date::<E, Saturn, _>(&jd, &ctx, E::saturn_heliocentric);
    let earth_opening = signed_opening(&pole_true, &planet_true);
    // Sunward unit vector from Saturn is opposite the heliocentric planet vector.
    let sun_opening = signed_opening_from_heliocentric(&pole_true, &planet_helio_true);
    SaturnRingOrientation {
        north_pole_position_angle: bearing_to_pole(&planet_true, &pole_true),
        ring_opening_earth: earth_opening,
        ring_opening_sun: sun_opening,
    }
}

impl Mars {
    /// Geocentric Mars physical orientation using [`DefaultEphemeris`].
    ///
    /// For an explicit backend call [`mars_orientation`].
    ///
    /// # Arguments
    ///
    /// * `jd` — observation epoch on the TT axis as a Julian Date.
    ///
    /// # Returns
    ///
    /// [`MarsOrientation`] for the default VSOP87 planetary backend.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::solar_system::Mars;
    ///
    /// let ori = Mars::orientation(siderust::J2000);
    /// assert!(ori.north_pole_position_angle.is_some());
    /// ```
    #[inline]
    pub fn orientation(jd: JulianDate) -> MarsOrientation {
        mars_orientation::<crate::coordinates::transform::context::DefaultEphemeris>(jd)
    }
}

impl Saturn {
    /// Geocentric Saturn ring orientation using [`DefaultEphemeris`].
    ///
    /// For an explicit backend call [`saturn_ring_orientation`].
    ///
    /// # Arguments
    ///
    /// * `jd` — observation epoch on the TT axis as a Julian Date.
    ///
    /// # Returns
    ///
    /// [`SaturnRingOrientation`] for the default VSOP87 planetary backend.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::solar_system::Saturn;
    ///
    /// let ori = Saturn::ring_orientation(siderust::J2000);
    /// assert!(ori.north_pole_position_angle.is_some());
    /// ```
    #[inline]
    pub fn ring_orientation(jd: JulianDate) -> SaturnRingOrientation {
        saturn_ring_orientation::<crate::coordinates::transform::context::DefaultEphemeris>(jd)
    }
}

type GeoEcliptic = cartesian::Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit>;
type GeoTrue = cartesian::Position<Geocentric, frames::EquatorialTrueOfDate, AstronomicalUnit>;
type HelioEcliptic = cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>;
type HelioTrue = cartesian::Position<Heliocentric, frames::EquatorialTrueOfDate, AstronomicalUnit>;

fn planet_sun_pole_true_of_date<E, Body, F>(
    jd: &JulianDate,
    ctx: &AstroContext<E>,
    planet_at: F,
) -> (
    GeoTrue,
    GeoTrue,
    spherical::Direction<frames::EquatorialTrueOfDate>,
    HelioTrue,
)
where
    E: PlanetaryEphemeris,
    Body: HasIauRotation,
    F: Fn(JulianDate) -> HelioEcliptic,
{
    let earth = E::earth_heliocentric(*jd);
    let (planet_geo, planet_helio, _jd_ret) = light_time_geocentric(jd, &earth, &planet_at);
    let sun_geo = GeoEcliptic::new(-earth.x(), -earth.y(), -earth.z());

    let planet_mean: cartesian::Position<
        Geocentric,
        frames::EquatorialMeanJ2000,
        AstronomicalUnit,
    > = TransformFrame::to_frame(&planet_geo);
    let sun_mean: cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, AstronomicalUnit> =
        TransformFrame::to_frame(&sun_geo);
    let helio_mean: cartesian::Position<
        Heliocentric,
        frames::EquatorialMeanJ2000,
        AstronomicalUnit,
    > = TransformFrame::to_frame(&planet_helio);

    let planet_true = planet_mean.to_frame_with::<frames::EquatorialTrueOfDate, _>(jd, ctx);
    let sun_true = sun_mean.to_frame_with::<frames::EquatorialTrueOfDate, _>(jd, ctx);
    let helio_true = helio_mean.to_frame_with::<frames::EquatorialTrueOfDate, _>(jd, ctx);
    let pole_true = iau_pole_true_of_date::<Body, E>(jd, ctx);

    (planet_true, sun_true, pole_true, helio_true)
}

fn light_time_geocentric<F>(
    jd: &JulianDate,
    earth: &HelioEcliptic,
    planet_at: &F,
) -> (GeoEcliptic, HelioEcliptic, JulianDate)
where
    F: Fn(JulianDate) -> HelioEcliptic,
{
    let mut jd_ret = *jd;
    let mut planet = planet_at(jd_ret);
    let mut geo = GeoEcliptic::new(
        planet.x() - earth.x(),
        planet.y() - earth.y(),
        planet.z() - earth.z(),
    );
    for _ in 0..LIGHT_TIME_ITERS {
        let delta_au = geo.distance().value();
        let tau_days = delta_au / AU_PER_DAY_C.value();
        jd_ret = JulianDate::new((jd.raw() - Days::new(tau_days)).value());
        planet = planet_at(jd_ret);
        geo = GeoEcliptic::new(
            planet.x() - earth.x(),
            planet.y() - earth.y(),
            planet.z() - earth.z(),
        );
    }
    (geo, planet, jd_ret)
}

fn iau_pole_true_of_date<Body, E>(
    jd: &JulianDate,
    ctx: &AstroContext<E>,
) -> spherical::Direction<frames::EquatorialTrueOfDate>
where
    Body: HasIauRotation,
    E: PlanetaryEphemeris,
{
    let params: &IauRotationParams = &Body::ROTATION;
    // IAU α₀, δ₀ are ICRF equatorial; EquatorialMeanJ2000 is the mean equator
    // of J2000 used throughout Siderust's planetary pipeline.
    let pole_j2000 = spherical::Direction::<frames::EquatorialMeanJ2000>::new(
        params.alpha0(*jd),
        params.delta0(*jd),
    );
    use crate::coordinates::transform::ext::SphericalDirectionAstroExt;
    pole_j2000.to_frame_with::<frames::EquatorialTrueOfDate, _>(jd, ctx)
}

fn bearing_to_position(from: &GeoTrue, toward: &GeoTrue) -> Option<Radians> {
    let from_dir = spherical::Position::from_cartesian(from).direction();
    let toward_dir = spherical::Position::from_cartesian(toward).direction();
    from_dir.position_angle_to(&toward_dir)
}

fn bearing_to_pole(
    from: &GeoTrue,
    pole: &spherical::Direction<frames::EquatorialTrueOfDate>,
) -> Option<Radians> {
    let from_dir = spherical::Position::from_cartesian(from).direction();
    from_dir.position_angle_to(pole)
}

fn signed_opening(
    pole: &spherical::Direction<frames::EquatorialTrueOfDate>,
    planet_geo: &GeoTrue,
) -> Radians {
    let Some(geo_dir) = planet_geo.direction() else {
        return Radians::new(0.0);
    };
    let pole_cart = pole.to_cartesian();
    // Planetocentric latitude of Earth: asin(n̂ · ê) with ê = planet → Earth = −ĝ.
    let sin_b = -(pole_cart.dot(&geo_dir));
    Radians::new(clamp_unit(sin_b).asin())
}

fn signed_opening_from_heliocentric(
    pole: &spherical::Direction<frames::EquatorialTrueOfDate>,
    planet_helio: &HelioTrue,
) -> Radians {
    let Some(helio_dir) = planet_helio.direction() else {
        return Radians::new(0.0);
    };
    let pole_cart = pole.to_cartesian();
    // Planetocentric latitude of the Sun: asin(n̂ · (−r̂_helio)).
    let sin_b = -(pole_cart.dot(&helio_dir));
    Radians::new(clamp_unit(sin_b).asin())
}

#[inline]
fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::{Barycentric, Geocentric, Heliocentric};
    use crate::coordinates::transform::context::DefaultEphemeris;
    use crate::ephemeris::{AuPerDay, Ephemeris, PlanetaryEphemeris};
    use crate::qtty::{AstronomicalUnit, Degree, Degrees, Kilometer};
    use core::sync::atomic::{AtomicUsize, Ordering};

    static MARS_CALLS: AtomicUsize = AtomicUsize::new(0);
    static EARTH_CALLS: AtomicUsize = AtomicUsize::new(0);

    /// Inject known true-of-date geometry through the ecliptic ephemeris contract.
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
        fn earth_heliocentric(
            _jd: JulianDate,
        ) -> cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>
        {
            EARTH_CALLS.fetch_add(1, Ordering::SeqCst);
            // Positions are anchored at J2000 so light-time epoch shifts do not
            // change the injected true-of-date geometry.
            let sun = source_position::<AstronomicalUnit>(crate::J2000, 90.0, 0.0);
            cartesian::Position::new(-sun.x(), -sun.y(), -sun.z())
        }

        fn moon_geocentric(
            _: JulianDate,
        ) -> cartesian::Position<Geocentric, frames::EclipticMeanJ2000, Kilometer> {
            panic!("unexpected Moon evaluation")
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

    impl<const CASE: u8> PlanetaryEphemeris for Cardinal<CASE> {
        fn mars_heliocentric(
            _jd: JulianDate,
        ) -> cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>
        {
            MARS_CALLS.fetch_add(1, Ordering::SeqCst);
            // Geocentric Mars = planet_helio − earth_helio. With Earth set so
            // −earth = Sun at (90°, 0°), choose planet so geo Mars is at (0°, 0°).
            // geo = planet − earth ⇒ planet = geo + earth = geo − sun_geo.
            let mars_geo = match CASE {
                0 | 1 => source_position::<AstronomicalUnit>(crate::J2000, 0.0, 0.0),
                _ => source_position::<AstronomicalUnit>(crate::J2000, 0.0, 90.0),
            };
            let sun_geo = source_position::<AstronomicalUnit>(crate::J2000, 90.0, 0.0);
            cartesian::Position::new(
                mars_geo.x() - sun_geo.x(),
                mars_geo.y() - sun_geo.y(),
                mars_geo.z() - sun_geo.z(),
            )
        }

        fn saturn_heliocentric(
            jd: JulianDate,
        ) -> cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>
        {
            Self::mars_heliocentric(jd)
        }
    }

    fn normalized(angle: Option<Radians>) -> Radians {
        let angle = angle.expect("defined bearing");
        assert!(angle.is_finite());
        assert!(angle >= Radians::new(0.0) && angle < Radians::FULL_TURN);
        angle
    }

    #[test]
    fn mars_bright_limb_cardinal_east_and_backend_calls() {
        // Mars at (0°,0°), Sun at (90°,0°) ⇒ bright limb PA = 90°.
        MARS_CALLS.store(0, Ordering::SeqCst);
        EARTH_CALLS.store(0, Ordering::SeqCst);
        let ori = mars_orientation::<Cardinal<0>>(crate::J2000);
        let bright = normalized(ori.bright_limb_position_angle);
        assert!((bright.to::<Degree>() - Degrees::new(90.0)).abs() < Degrees::new(1e-8));
        // Light-time uses Earth once + planet LIGHT_TIME_ITERS+1 times (initial + iters).
        assert_eq!(EARTH_CALLS.load(Ordering::SeqCst), 1);
        assert_eq!(
            MARS_CALLS.load(Ordering::SeqCst),
            usize::from(LIGHT_TIME_ITERS) + 1
        );
    }

    #[test]
    fn convenience_wrappers_match_default_backend() {
        let jd = crate::J2000;
        assert_eq!(
            Mars::orientation(jd),
            mars_orientation::<DefaultEphemeris>(jd)
        );
        assert_eq!(
            Saturn::ring_orientation(jd),
            saturn_ring_orientation::<DefaultEphemeris>(jd)
        );
    }

    #[test]
    fn meeus_42a_mars_position_angles() {
        // Meeus, Astronomical Algorithms 2nd ed., Example 42.a, p. 291.
        // Independently reproduced by soniakeys/meeus (Physical):
        // JDE 2448935.500683 → P = 347.64°, Q (defect) = 279.91°.
        // Bright limb = Q + 180° → 99.91°.
        // Tolerances absorb IAU pole-model and VSOP87/FK5 differences vs Meeus.
        let jd = JulianDate::new(2_448_935.500_683);
        let ori = mars_orientation::<DefaultEphemeris>(jd);
        let p = normalized(ori.north_pole_position_angle).to::<Degree>();
        let bright = normalized(ori.bright_limb_position_angle).to::<Degree>();
        let defect = (ori.bright_limb_position_angle.unwrap() + Radians::HALF_TURN)
            .normalize()
            .to::<Degree>();
        // Bright-limb/defect track Meeus closely; pole PA differs by <1° because
        // Siderust uses the modern IAU linear pole vs Meeus's ecliptic λ₀/β₀ model.
        assert!(
            (p - Degrees::new(347.64)).abs() < Degrees::new(1.0),
            "pole PA {p:?} vs Meeus 347.64°"
        );
        assert!(
            (bright - Degrees::new(99.91)).abs() < Degrees::new(0.05),
            "bright PA {bright:?} vs 99.91°"
        );
        assert!(
            (defect - Degrees::new(279.91)).abs() < Degrees::new(0.05),
            "defect PA {defect:?} vs Meeus Q 279.91°"
        );
    }

    #[test]
    fn meeus_45a_saturn_rings() {
        // Meeus Example 45.a, p. 320 / soniakeys/meeus saturnring:
        // JDE 2448972.50068 → B = 16.442°, B′ = 14.679°, P = 6.741°.
        let jd = JulianDate::new(2_448_972.500_68);
        let ori = saturn_ring_orientation::<DefaultEphemeris>(jd);
        let p = normalized(ori.north_pole_position_angle).to::<Degree>();
        let b = ori.ring_opening_earth.to::<Degree>();
        let bp = ori.ring_opening_sun.to::<Degree>();
        assert!(
            (p - Degrees::new(6.741)).abs() < Degrees::new(0.05),
            "pole/ring PA {p:?} vs Meeus 6.741°"
        );
        assert!(
            (b - Degrees::new(16.442)).abs() < Degrees::new(0.05),
            "B {b:?} vs Meeus 16.442°"
        );
        assert!(
            (bp - Degrees::new(14.679)).abs() < Degrees::new(0.05),
            "B′ {bp:?} vs Meeus 14.679°"
        );
        assert!(b.value() > 0.0 && bp.value() > 0.0);
    }

    #[test]
    fn saturn_opening_sign_and_bounds() {
        for jd in [
            crate::J2000,
            JulianDate::new(2_448_972.500_68),
            JulianDate::new(2_460_000.0),
        ] {
            let ori = Saturn::ring_orientation(jd);
            let half_pi = Radians::new(core::f64::consts::FRAC_PI_2);
            assert!(ori.ring_opening_earth.abs() <= half_pi + Radians::new(1e-12));
            assert!(ori.ring_opening_sun.abs() <= half_pi + Radians::new(1e-12));
        }
    }

    #[test]
    fn angle_normalization_and_undefined_geometry() {
        let origin = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
        );
        assert_eq!(
            origin.position_angle_to(&origin),
            None,
            "coincident directions are undefined"
        );
        let pole = spherical::Direction::<frames::EquatorialTrueOfDate>::new(
            Degrees::new(0.0),
            Degrees::new(90.0),
        );
        assert_eq!(pole.position_angle_to(&origin), None);
    }

    #[test]
    fn multi_epoch_finite_vsop87() {
        for jd in [
            JulianDate::new(2_448_935.500_683),
            JulianDate::new(2_448_972.500_68),
            crate::J2000,
            JulianDate::new(2_460_000.5),
        ] {
            let mars = Mars::orientation(jd);
            assert!(mars.north_pole_position_angle.unwrap().is_finite());
            assert!(mars.bright_limb_position_angle.unwrap().is_finite());
            let sat = Saturn::ring_orientation(jd);
            assert!(sat.north_pole_position_angle.unwrap().is_finite());
            assert!(sat.ring_opening_earth.is_finite());
            assert!(sat.ring_opening_sun.is_finite());
        }
    }
}
