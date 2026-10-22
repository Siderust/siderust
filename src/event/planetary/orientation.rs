// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Planetary physical orientation
//!
//! ## Scientific scope
//!
//! Light-time geometric sky-plane orientation of a major planet with an IAU
//! WGCCRE rotational model:
//!
//! - Position angle of the projected north rotational pole
//! - Bright-limb position angle (Sunward bearing)
//! - For Saturn, signed ring opening angles from Earth and from the Sun
//!
//! Position angles are measured from celestial north through east in a
//! selected equatorial [`SkyOrientationFrame`], normalised to `[0, 2π)` when
//! defined. Signed ring openings are saturnicentric latitudes of the Earth and
//! Sun relative to the ring plane (IAU north pole); **positive** means the
//! north face is presented to that observer. Openings are frame-invariant.
//!
//! The IAU pole (`α₀`, `δ₀`) is interpreted in **ICRS/ICRF** and transformed
//! into the sky frame. Planetary heliocentric states are evaluated at a
//! light-time retarded TT epoch (two fixed-point iterations); Earth is taken
//! at the observation epoch. Both share backend `E: PlanetState<Body>`.
//!
//! These are **not** fully apparent places: annual/diurnal aberration,
//! gravitational deflection, and atmospheric refraction are omitted. The
//! topocentric path applies parallax to the planet and Sun with a shared site.
//! TT is the working scale; TT−TDB is neglected for the linear pole rates.
//!
//! ## Technical scope
//!
//! - [`body_orientation`] — generic over body, ephemeris, observation center,
//!   and sky frame
//! - [`PlanetaryObservationCenter`] — [`Geocentric`] / [`Topocentric`] inputs
//! - [`SkyOrientationFrame`] — [`EquatorialTrueOfDate`] (default),
//!   [`EquatorialMeanJ2000`], [`ICRS`]
//! - Convenience wrappers [`mars_orientation`], [`saturn_ring_orientation`],
//!   [`Mars::orientation`], [`Saturn::ring_orientation`]
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

use crate::astro::HasIauRotation;
use crate::bodies::solar_system::{Mars, Saturn};
use crate::coordinates::centers::{
    Geocentric, Geodetic, Heliocentric, ReferenceCenter, Topocentric,
};
use crate::coordinates::frames::{
    self, EquatorialMeanJ2000, EquatorialTrueOfDate, MutableFrame, ReferenceFrame, ECEF, ICRS,
};
use crate::coordinates::transform::centers::position::to_topocentric::to_topocentric_with;
use crate::coordinates::transform::ext::SphericalDirectionAstroExt;
use crate::coordinates::transform::providers::frame_rotation_with;
use crate::coordinates::transform::{AstroContext, TransformContext, TransformFrame};
use crate::coordinates::{cartesian, spherical};
use crate::ephemeris::{Ephemeris, PlanetState};
use crate::qtty::{velocity::AU_PER_DAY_C, AstronomicalUnit, Days, LengthUnit, Radians};
use crate::time::JulianDate;

#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};

/// Fixed light-time iterations (Meeus planetary physical ephemerides use two).
const LIGHT_TIME_ITERS: u8 = 2;

// ---------------------------------------------------------------------------
// Sky frame
// ---------------------------------------------------------------------------

/// Equatorial frames that admit a celestial north-through-east position-angle
/// convention for disc orientation.
///
/// Not every [`ReferenceFrame`] is a sky-orientation frame. Implementations
/// convert ICRS pole directions and mean-equatorial J2000 positions into `Self`.
pub trait SkyOrientationFrame: MutableFrame + Sized {
    /// Rotate an ICRS direction into this sky frame.
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self>;

    /// Rotate a mean-equatorial J2000 position into this sky frame.
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy;

    /// Rotate a topocentric mean-equatorial J2000 position into this sky frame.
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>;
}

impl SkyOrientationFrame for EquatorialTrueOfDate {
    #[inline]
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        dir.to_frame_with(jd, ctx)
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy,
    {
        rotate_mean_equatorial_position(pos, jd, ctx)
    }

    #[inline]
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>,
    {
        let rotation = frame_rotation_with::<EquatorialMeanJ2000, Self, _>(*jd, ctx);
        cartesian::Position::from_array(*pos.center_params(), rotation * *pos.as_array())
    }
}

impl SkyOrientationFrame for EquatorialMeanJ2000 {
    #[inline]
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        dir.to_frame_with(jd, ctx)
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        _jd: &JulianDate,
        _ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy,
    {
        TransformFrame::to_frame(pos)
    }

    #[inline]
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        _jd: &JulianDate,
        _ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>,
    {
        TransformFrame::to_frame(pos)
    }
}

impl SkyOrientationFrame for ICRS {
    #[inline]
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        _jd: &JulianDate,
        _ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        *dir
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy,
    {
        rotate_mean_equatorial_position(pos, jd, ctx)
    }

    #[inline]
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>,
    {
        let rotation = frame_rotation_with::<EquatorialMeanJ2000, Self, _>(*jd, ctx);
        cartesian::Position::from_array(*pos.center_params(), rotation * *pos.as_array())
    }
}

/// Rotate a mean-equatorial J2000 position into `F` without requiring
/// `ReferenceCenter::Params = ()` (unlike `PositionAstroExt`).
fn rotate_mean_equatorial_position<C, U, Eph, F>(
    pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
    jd: &JulianDate,
    ctx: &AstroContext<Eph>,
) -> cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    U: LengthUnit,
    Eph: Ephemeris,
    F: MutableFrame,
    C::Params: Copy,
    (): crate::coordinates::transform::providers::FrameRotationProvider<EquatorialMeanJ2000, F>,
{
    let rotation = frame_rotation_with::<EquatorialMeanJ2000, F, _>(*jd, ctx);
    cartesian::Position::from_array(*pos.center_params(), rotation * *pos.as_array())
}

// ---------------------------------------------------------------------------
// Observation center
// ---------------------------------------------------------------------------

/// Observation origin for planetary physical orientation.
///
/// Implemented for [`Geocentric`] and [`Topocentric`]. Dispatch is
/// monomorphized; unsupported centers fail to compile.
pub trait PlanetaryObservationCenter: ReferenceCenter {
    /// Observation input: TT [`JulianDate`] for [`Geocentric`], or
    /// `(JulianDate, Geodetic<ECEF>)` for [`Topocentric`].
    type Input;

    /// Body orientation in sky frame `F` with backend `E`.
    fn body_orientation<Body, E, F>(input: Self::Input) -> BodyOrientation
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame;

    /// Saturn ring orientation in sky frame `F` with backend `E`.
    fn saturn_ring_orientation<E, F>(input: Self::Input) -> SaturnRingOrientation
    where
        E: PlanetState<Saturn>,
        F: SkyOrientationFrame;
}

impl PlanetaryObservationCenter for Geocentric {
    type Input = JulianDate;

    fn body_orientation<Body, E, F>(jd: Self::Input) -> BodyOrientation
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let geo = light_time_states::<E, Body>(&jd);
        let mean = to_mean_equatorial_geocentric(&geo);
        let sky = SkyBundle {
            planet: F::position_from_mean_equatorial(&mean.planet, &jd, &ctx),
            sun: F::position_from_mean_equatorial(&mean.sun, &jd, &ctx),
            pole: iau_pole_in_frame::<Body, E, F>(&jd, &ctx),
            planet_helio: F::position_from_mean_equatorial(&mean.planet_helio, &jd, &ctx),
        };
        body_from_sky(&sky)
    }

    fn saturn_ring_orientation<E, F>(jd: Self::Input) -> SaturnRingOrientation
    where
        E: PlanetState<Saturn>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let geo = light_time_states::<E, Saturn>(&jd);
        let mean = to_mean_equatorial_geocentric(&geo);
        let sky = SkyBundle {
            planet: F::position_from_mean_equatorial(&mean.planet, &jd, &ctx),
            sun: F::position_from_mean_equatorial(&mean.sun, &jd, &ctx),
            pole: iau_pole_in_frame::<Saturn, E, F>(&jd, &ctx),
            planet_helio: F::position_from_mean_equatorial(&mean.planet_helio, &jd, &ctx),
        };
        saturn_from_sky(&sky)
    }
}

impl PlanetaryObservationCenter for Topocentric {
    type Input = (JulianDate, Geodetic<ECEF>);

    fn body_orientation<Body, E, F>((jd, site): Self::Input) -> BodyOrientation
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let geo = light_time_states::<E, Body>(&jd);
        let mean = to_mean_equatorial_geocentric(&geo);
        let planet_topo = to_topocentric_with(&mean.planet, site, jd, &ctx);
        let sun_topo = to_topocentric_with(&mean.sun, site, jd, &ctx);
        let sky = SkyBundle {
            planet: F::topocentric_from_mean_equatorial::<E, _>(&planet_topo, &jd, &ctx),
            sun: F::topocentric_from_mean_equatorial::<E, _>(&sun_topo, &jd, &ctx),
            pole: iau_pole_in_frame::<Body, E, F>(&jd, &ctx),
            planet_helio: F::position_from_mean_equatorial(&mean.planet_helio, &jd, &ctx),
        };
        body_from_sky(&sky)
    }

    fn saturn_ring_orientation<E, F>((jd, site): Self::Input) -> SaturnRingOrientation
    where
        E: PlanetState<Saturn>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let geo = light_time_states::<E, Saturn>(&jd);
        let mean = to_mean_equatorial_geocentric(&geo);
        let planet_topo = to_topocentric_with(&mean.planet, site, jd, &ctx);
        let sun_topo = to_topocentric_with(&mean.sun, site, jd, &ctx);
        let sky = SkyBundle {
            planet: F::topocentric_from_mean_equatorial::<E, _>(&planet_topo, &jd, &ctx),
            sun: F::topocentric_from_mean_equatorial::<E, _>(&sun_topo, &jd, &ctx),
            pole: iau_pole_in_frame::<Saturn, E, F>(&jd, &ctx),
            planet_helio: F::position_from_mean_equatorial(&mean.planet_helio, &jd, &ctx),
        };
        saturn_from_sky(&sky)
    }
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Physical orientation of a planetary disc on the sky.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyOrientation {
    /// Position angle of the north rotational pole, north through east in
    /// `[0, 2π)`, or `None` when affn reports undefined sky-plane geometry.
    pub north_pole_position_angle: Option<Radians>,
    /// Bright-limb position angle (bearing toward the Sun), same convention, or
    /// `None` when undefined. Greatest defect of illumination is
    /// `(bright + Radians::HALF_TURN).normalize()` when `bright` is defined.
    pub bright_limb_position_angle: Option<Radians>,
}

/// Mars-specific alias of [`BodyOrientation`] (stable name for callers).
pub type MarsOrientation = BodyOrientation;

/// Saturn ring / pole orientation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SaturnRingOrientation {
    /// Position angle of the projected north pole (ring axis), or `None` when
    /// undefined.
    pub north_pole_position_angle: Option<Radians>,
    /// Signed ring opening as seen from the observation center (saturnicentric
    /// latitude of the observer). Positive ⇒ north face toward the observer.
    pub ring_opening_earth: Radians,
    /// Signed ring opening as seen from the Sun. Positive ⇒ north face
    /// illuminated. Independent of observation center and sky frame.
    pub ring_opening_sun: Radians,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Planetary disc orientation with compile-time body, backend, center, and frame.
///
/// # Arguments
///
/// * `input` — TT [`JulianDate`] for [`Geocentric`], or
///   `(JulianDate, Geodetic<ECEF>)` for [`Topocentric`].
///
/// # Returns
///
/// [`BodyOrientation`]. See the [module documentation](self) for corrections.
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::solar_system::Mars;
/// use siderust::coordinates::centers::Geocentric;
/// use siderust::coordinates::frames::EquatorialTrueOfDate;
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::body_orientation;
///
/// let ori = body_orientation::<Mars, Vsop87Ephemeris, Geocentric, EquatorialTrueOfDate>(
///     siderust::J2000,
/// );
/// assert!(ori.north_pole_position_angle.unwrap().is_finite());
/// ```
///
/// Unsupported centers cannot be selected:
///
/// ```compile_fail,E0277
/// use siderust::bodies::solar_system::Mars;
/// use siderust::coordinates::centers::Heliocentric;
/// use siderust::coordinates::frames::EquatorialTrueOfDate;
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::body_orientation;
/// let _ = body_orientation::<Mars, Vsop87Ephemeris, Heliocentric, EquatorialTrueOfDate>(
///     siderust::J2000,
/// );
/// ```
pub fn body_orientation<Body, E, C, F>(input: C::Input) -> BodyOrientation
where
    Body: HasIauRotation,
    E: PlanetState<Body>,
    C: PlanetaryObservationCenter,
    F: SkyOrientationFrame,
{
    C::body_orientation::<Body, E, F>(input)
}

/// Geocentric Mars orientation in the true equator/equinox of date.
///
/// Equivalent to
/// `body_orientation::<Mars, E, Geocentric, EquatorialTrueOfDate>(jd)`.
///
/// # Examples
///
/// ```rust
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::mars_orientation;
///
/// let ori = mars_orientation::<Vsop87Ephemeris>(siderust::J2000);
/// assert!(ori.north_pole_position_angle.unwrap().is_finite());
/// ```
#[inline]
pub fn mars_orientation<E: PlanetState<Mars>>(jd: JulianDate) -> MarsOrientation {
    body_orientation::<Mars, E, Geocentric, EquatorialTrueOfDate>(jd)
}

/// Saturn ring orientation with compile-time backend, center, and sky frame.
///
/// # Examples
///
/// ```rust
/// use siderust::coordinates::centers::Geocentric;
/// use siderust::coordinates::frames::EquatorialTrueOfDate;
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::saturn_rings;
/// use siderust::qtty::Degrees;
///
/// let ori = saturn_rings::<Vsop87Ephemeris, Geocentric, EquatorialTrueOfDate>(siderust::J2000);
/// assert!(ori.ring_opening_earth.to::<siderust::qtty::Degree>().abs() < Degrees::new(90.0));
/// ```
pub fn saturn_rings<E, C, F>(input: C::Input) -> SaturnRingOrientation
where
    E: PlanetState<Saturn>,
    C: PlanetaryObservationCenter,
    F: SkyOrientationFrame,
{
    C::saturn_ring_orientation::<E, F>(input)
}

/// Geocentric Saturn ring orientation in the true equator/equinox of date.
///
/// Equivalent to
/// `saturn_rings::<E, Geocentric, EquatorialTrueOfDate>(jd)`.
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
#[inline]
pub fn saturn_ring_orientation<E: PlanetState<Saturn>>(jd: JulianDate) -> SaturnRingOrientation {
    saturn_rings::<E, Geocentric, EquatorialTrueOfDate>(jd)
}

impl Mars {
    /// Geocentric Mars orientation using [`DefaultEphemeris`] in
    /// [`EquatorialTrueOfDate`].
    ///
    /// For explicit center/frame/backend selection see [`orientation_at`](Self::orientation_at)
    /// and [`body_orientation`].
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
        Self::orientation_at::<Geocentric, EquatorialTrueOfDate>(jd)
    }

    /// Mars orientation with default ephemeris and explicit center/frame.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::solar_system::Mars;
    /// use siderust::coordinates::centers::Geocentric;
    /// use siderust::coordinates::frames::ICRS;
    ///
    /// let ori = Mars::orientation_at::<Geocentric, ICRS>(siderust::J2000);
    /// assert!(ori.north_pole_position_angle.is_some());
    /// ```
    #[inline]
    pub fn orientation_at<C, F>(input: C::Input) -> MarsOrientation
    where
        C: PlanetaryObservationCenter,
        F: SkyOrientationFrame,
    {
        body_orientation::<Mars, crate::coordinates::transform::context::DefaultEphemeris, C, F>(
            input,
        )
    }
}

impl Saturn {
    /// Geocentric Saturn ring orientation using [`DefaultEphemeris`] in
    /// [`EquatorialTrueOfDate`].
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
        Self::ring_orientation_at::<Geocentric, EquatorialTrueOfDate>(jd)
    }

    /// Saturn ring orientation with default ephemeris and explicit center/frame.
    #[inline]
    pub fn ring_orientation_at<C, F>(input: C::Input) -> SaturnRingOrientation
    where
        C: PlanetaryObservationCenter,
        F: SkyOrientationFrame,
    {
        saturn_rings::<crate::coordinates::transform::context::DefaultEphemeris, C, F>(input)
    }
}

// ---------------------------------------------------------------------------
// Internal geometry
// ---------------------------------------------------------------------------

type HelioEcliptic = cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>;
type GeoEcliptic = cartesian::Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit>;
type GeoMean = cartesian::Position<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>;
type HelioMean = cartesian::Position<Heliocentric, EquatorialMeanJ2000, AstronomicalUnit>;

struct LightTimeStates {
    planet: GeoEcliptic,
    sun: GeoEcliptic,
    planet_helio: HelioEcliptic,
}

struct MeanEquatorial {
    planet: GeoMean,
    sun: GeoMean,
    planet_helio: HelioMean,
}

struct SkyBundle<C: ReferenceCenter, F: ReferenceFrame> {
    planet: cartesian::Position<C, F, AstronomicalUnit>,
    sun: cartesian::Position<C, F, AstronomicalUnit>,
    pole: spherical::Direction<F>,
    planet_helio: cartesian::Position<Heliocentric, F, AstronomicalUnit>,
}

fn light_time_states<E, Body>(jd: &JulianDate) -> LightTimeStates
where
    E: PlanetState<Body>,
    Body: HasIauRotation,
{
    let earth = E::earth_heliocentric(*jd);
    let (planet, planet_helio, _) = light_time_geocentric(jd, &earth, &|t| E::heliocentric(t));
    let sun = GeoEcliptic::new(-earth.x(), -earth.y(), -earth.z());
    LightTimeStates {
        planet,
        sun,
        planet_helio,
    }
}

fn to_mean_equatorial_geocentric(geo: &LightTimeStates) -> MeanEquatorial {
    MeanEquatorial {
        planet: TransformFrame::to_frame(&geo.planet),
        sun: TransformFrame::to_frame(&geo.sun),
        planet_helio: TransformFrame::to_frame(&geo.planet_helio),
    }
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

fn iau_pole_in_frame<Body, E, F>(jd: &JulianDate, ctx: &AstroContext<E>) -> spherical::Direction<F>
where
    Body: HasIauRotation,
    E: Ephemeris,
    F: SkyOrientationFrame,
{
    let params = &Body::ROTATION;
    let pole_icrs = spherical::Direction::<ICRS>::new(params.alpha0(*jd), params.delta0(*jd));
    F::direction_from_icrs(&pole_icrs, jd, ctx)
}

fn body_from_sky<C, F>(sky: &SkyBundle<C, F>) -> BodyOrientation
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    BodyOrientation {
        north_pole_position_angle: bearing_to_pole(&sky.planet, &sky.pole),
        bright_limb_position_angle: bearing_to_position(&sky.planet, &sky.sun),
    }
}

fn saturn_from_sky<C, F>(sky: &SkyBundle<C, F>) -> SaturnRingOrientation
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    SaturnRingOrientation {
        north_pole_position_angle: bearing_to_pole(&sky.planet, &sky.pole),
        ring_opening_earth: signed_opening(&sky.pole, &sky.planet),
        ring_opening_sun: signed_opening(&sky.pole, &sky.planet_helio),
    }
}

fn bearing_to_position<C1, C2, F>(
    from: &cartesian::Position<C1, F, AstronomicalUnit>,
    toward: &cartesian::Position<C2, F, AstronomicalUnit>,
) -> Option<Radians>
where
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F: ReferenceFrame,
{
    let from_dir = spherical::Position::from_cartesian(from).direction();
    let toward_dir = spherical::Position::from_cartesian(toward).direction();
    from_dir.position_angle_to(&toward_dir)
}

fn bearing_to_pole<C, F>(
    from: &cartesian::Position<C, F, AstronomicalUnit>,
    pole: &spherical::Direction<F>,
) -> Option<Radians>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let from_dir = spherical::Position::from_cartesian(from).direction();
    from_dir.position_angle_to(pole)
}

fn signed_opening<C, F>(
    pole: &spherical::Direction<F>,
    planet: &cartesian::Position<C, F, AstronomicalUnit>,
) -> Radians
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let Some(dir) = planet.direction() else {
        return Radians::new(0.0);
    };
    let pole_cart = pole.to_cartesian();
    // Planetocentric latitude of the observer: asin(n̂ · (−r̂)).
    let sin_b = -(pole_cart.dot(&dir));
    Radians::new(clamp_unit(sin_b).asin())
}

#[inline]
fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Jupiter;
    use crate::coordinates::centers::{Barycentric, Geocentric, Heliocentric};
    use crate::coordinates::transform::context::DefaultEphemeris;
    use crate::coordinates::transform::ext::PositionAstroExt;
    use crate::ephemeris::{AuPerDay, Ephemeris, PlanetState};
    use crate::qtty::{AstronomicalUnit, Degree, Degrees, Kilometer, Meters};
    use core::sync::atomic::{AtomicUsize, Ordering};

    static MARS_CALLS: AtomicUsize = AtomicUsize::new(0);
    static EARTH_CALLS: AtomicUsize = AtomicUsize::new(0);

    struct Cardinal<const CASE: u8>;

    fn source_position<U: crate::qtty::LengthUnit>(
        jd: JulianDate,
        ra: f64,
        dec: f64,
    ) -> cartesian::Position<Geocentric, frames::EclipticMeanJ2000, U> {
        let true_pos = spherical::Position::<Geocentric, EquatorialTrueOfDate, U>::new(
            Degrees::new(ra),
            Degrees::new(dec),
            crate::qtty::Quantity::<U>::new(1.0),
        )
        .to_cartesian();
        let mean = true_pos.to_frame_with::<EquatorialMeanJ2000, _>(
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

    impl<const CASE: u8> PlanetState<Mars> for Cardinal<CASE> {
        fn heliocentric(
            _jd: JulianDate,
        ) -> cartesian::Position<Heliocentric, frames::EclipticMeanJ2000, AstronomicalUnit>
        {
            MARS_CALLS.fetch_add(1, Ordering::SeqCst);
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
    }

    fn normalized(angle: Option<Radians>) -> Radians {
        let angle = angle.expect("defined bearing");
        assert!(angle.is_finite());
        assert!(angle >= Radians::new(0.0) && angle < Radians::FULL_TURN);
        angle
    }

    #[test]
    fn mars_bright_limb_cardinal_east_and_backend_calls() {
        MARS_CALLS.store(0, Ordering::SeqCst);
        EARTH_CALLS.store(0, Ordering::SeqCst);
        let ori = mars_orientation::<Cardinal<0>>(crate::J2000);
        let bright = normalized(ori.bright_limb_position_angle);
        assert!((bright.to::<Degree>() - Degrees::new(90.0)).abs() < Degrees::new(1e-8));
        assert_eq!(EARTH_CALLS.load(Ordering::SeqCst), 1);
        assert_eq!(
            MARS_CALLS.load(Ordering::SeqCst),
            usize::from(LIGHT_TIME_ITERS) + 1
        );
    }

    #[test]
    fn convenience_wrappers_match_generic_defaults() {
        let jd = crate::J2000;
        assert_eq!(
            Mars::orientation(jd),
            mars_orientation::<DefaultEphemeris>(jd)
        );
        assert_eq!(
            Mars::orientation(jd),
            body_orientation::<Mars, DefaultEphemeris, Geocentric, EquatorialTrueOfDate>(jd)
        );
        assert_eq!(
            Saturn::ring_orientation(jd),
            saturn_ring_orientation::<DefaultEphemeris>(jd)
        );
        assert_eq!(
            Saturn::ring_orientation(jd),
            saturn_rings::<DefaultEphemeris, Geocentric, EquatorialTrueOfDate>(jd)
        );
    }

    #[test]
    fn meeus_42a_mars_position_angles() {
        let jd = JulianDate::new(2_448_935.500_683);
        let ori = mars_orientation::<DefaultEphemeris>(jd);
        let p = normalized(ori.north_pole_position_angle).to::<Degree>();
        let bright = normalized(ori.bright_limb_position_angle).to::<Degree>();
        let defect = (ori.bright_limb_position_angle.unwrap() + Radians::HALF_TURN)
            .normalize()
            .to::<Degree>();
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
    fn jupiter_generic_body_is_finite() {
        let ori = body_orientation::<Jupiter, DefaultEphemeris, Geocentric, EquatorialTrueOfDate>(
            crate::J2000,
        );
        assert!(ori.north_pole_position_angle.unwrap().is_finite());
        assert!(ori.bright_limb_position_angle.unwrap().is_finite());
    }

    #[test]
    fn geocentric_vs_topocentric_differ_and_site_matters() {
        let jd = JulianDate::new(2_448_935.500_683);
        let geo = Mars::orientation(jd);
        let sites = [
            Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0)),
            Geodetic::<ECEF>::new(Degrees::new(-70.4), Degrees::new(-24.6), Meters::new(0.0)),
        ];
        let mut topocentric = [Radians::new(0.0); 2];
        for (i, site) in sites.iter().enumerate() {
            let topo = Mars::orientation_at::<Topocentric, EquatorialTrueOfDate>((jd, *site));
            assert!(
                (geo.north_pole_position_angle.unwrap() - topo.north_pole_position_angle.unwrap())
                    .abs()
                    > Radians::new(1e-8)
            );
            topocentric[i] = topo.north_pole_position_angle.unwrap();
        }
        assert!((topocentric[0] - topocentric[1]).abs() > Radians::new(1e-8));
    }

    #[test]
    fn equatorial_frames_transform_and_openings_are_invariant() {
        let jd = JulianDate::new(2_448_972.500_68);
        let tod = saturn_rings::<DefaultEphemeris, Geocentric, EquatorialTrueOfDate>(jd);
        let mean = saturn_rings::<DefaultEphemeris, Geocentric, EquatorialMeanJ2000>(jd);
        let icrs = saturn_rings::<DefaultEphemeris, Geocentric, ICRS>(jd);
        // Ring openings are physical, not sky-plane, and must match across frames.
        assert!((tod.ring_opening_earth - mean.ring_opening_earth).abs() < Radians::new(1e-10));
        assert!((tod.ring_opening_earth - icrs.ring_opening_earth).abs() < Radians::new(1e-10));
        assert!((tod.ring_opening_sun - mean.ring_opening_sun).abs() < Radians::new(1e-10));
        assert!((tod.ring_opening_sun - icrs.ring_opening_sun).abs() < Radians::new(1e-10));
        // Sky-plane PAs generally differ between true-of-date and mean/ICRS.
        assert!(
            (tod.north_pole_position_angle.unwrap() - mean.north_pole_position_angle.unwrap())
                .abs()
                > Radians::new(1e-6)
                || (tod.north_pole_position_angle.unwrap()
                    - icrs.north_pole_position_angle.unwrap())
                .abs()
                    > Radians::new(1e-6)
        );
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
        let origin =
            spherical::Direction::<EquatorialTrueOfDate>::new(Degrees::new(0.0), Degrees::new(0.0));
        assert_eq!(origin.position_angle_to(&origin), None);
        let pole = spherical::Direction::<EquatorialTrueOfDate>::new(
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
