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
//! - For Saturn, signed ring opening angles from the observer and from the Sun
//!
//! Position angles are measured from celestial north through east in a
//! selected equatorial [`SkyOrientationFrame`], normalised to `[0, 2π)` when
//! defined. Signed ring openings are saturnicentric latitudes of the observer
//! and Sun relative to the ring plane (IAU north pole); **positive** means the
//! north face is presented to that body. Openings are frame-invariant.
//!
//! ### Epochs and corrections
//!
//! | Quantity | Epoch |
//! |----------|-------|
//! | Observer / Earth heliocentric state | Reception (TT observation epoch) |
//! | Planet heliocentric state | Emission (light-time retarded TT) |
//! | IAU pole `α₀`, `δ₀` | Emission (retarded TT; rates are TDB-based, TT−TDB neglected) |
//! | Frame rotations (precession/nutation) | Reception |
//!
//! Light-time uses the **observer–planet** range (geocentric distance, or
//! topocentric distance including the site vector), iterated twice. The Sun
//! for bright-limb / solar opening is the geometric geocentric Sun at
//! reception (`−r_Earth`), consistent with a single shared ephemeris backend.
//!
//! Applied corrections match [`CorrectionPolicy::LIGHT_TIME`] only: retarded
//! planetary geometry and topocentric parallax. **Omitted:** annual/diurnal
//! aberration, gravitational deflection, atmospheric refraction. Results are
//! therefore **light-time geometric**, not fully apparent places. Do not label
//! them as apparent without those corrections. Near conjunction or celestial
//! poles, sky-plane bearings may be `None`; a vanishing planet vector yields
//! `None` openings (zero opening remains a valid `Some(0)` edge-on case).
//!
//! ## Technical scope
//!
//! - [`body_orientation`] — generic over body, ephemeris, observation center,
//!   and sky frame
//! - [`PlanetaryObservationCenter`] — [`Geocentric`] / [`Topocentric`] inputs
//! - [`SkyOrientationFrame`] — [`EquatorialTrueOfDate`] (default),
//!   [`EquatorialMeanJ2000`], [`ICRS`]
//! - Shared observation-geometry pipeline; Saturn rings layer openings on top
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

use crate::astro::apparent::CorrectionPolicy;
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

/// Correction stages applied by this module (light-time only).
pub const ORIENTATION_CORRECTIONS: CorrectionPolicy = CorrectionPolicy::LIGHT_TIME;

// ---------------------------------------------------------------------------
// Sky frame
// ---------------------------------------------------------------------------

/// Equatorial frames that admit a celestial north-through-east position-angle
/// convention for disc orientation.
///
/// Not every [`ReferenceFrame`] is a sky-orientation frame. Implementations
/// convert ICRS pole directions and mean-equatorial J2000 positions into `Self`.
pub trait SkyOrientationFrame: MutableFrame + Sized {
    /// Rotate an ICRS direction into this sky frame at the reception epoch.
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        jd_reception: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self>;

    /// Rotate a mean-equatorial J2000 position into this sky frame.
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd_reception: &JulianDate,
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
        jd_reception: &JulianDate,
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
        jd_reception: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        dir.to_frame_with(jd_reception, ctx)
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd_reception: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy,
    {
        rotate_mean_equatorial_position(pos, jd_reception, ctx)
    }

    #[inline]
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        jd_reception: &JulianDate,
        ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>,
    {
        let rotation = frame_rotation_with::<EquatorialMeanJ2000, Self, _>(*jd_reception, ctx);
        cartesian::Position::from_array(*pos.center_params(), rotation * *pos.as_array())
    }
}

impl SkyOrientationFrame for EquatorialMeanJ2000 {
    #[inline]
    fn direction_from_icrs<Eph: Ephemeris>(
        dir: &spherical::Direction<ICRS>,
        jd_reception: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        dir.to_frame_with(jd_reception, ctx)
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        _jd_reception: &JulianDate,
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
        _jd_reception: &JulianDate,
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
        _jd_reception: &JulianDate,
        _ctx: &AstroContext<Eph>,
    ) -> spherical::Direction<Self> {
        *dir
    }

    #[inline]
    fn position_from_mean_equatorial<C, U, Eph>(
        pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
        jd_reception: &JulianDate,
        ctx: &AstroContext<Eph>,
    ) -> cartesian::Position<C, Self, U>
    where
        C: ReferenceCenter,
        U: LengthUnit,
        Eph: Ephemeris,
        C::Params: Copy,
    {
        rotate_mean_equatorial_position(pos, jd_reception, ctx)
    }

    #[inline]
    fn topocentric_from_mean_equatorial<Eph, Ctx>(
        pos: &cartesian::Position<Topocentric, EquatorialMeanJ2000, AstronomicalUnit>,
        jd_reception: &JulianDate,
        ctx: &Ctx,
    ) -> cartesian::Position<Topocentric, Self, AstronomicalUnit>
    where
        Eph: Ephemeris,
        Ctx: TransformContext<Eph = Eph>,
    {
        let rotation = frame_rotation_with::<EquatorialMeanJ2000, Self, _>(*jd_reception, ctx);
        cartesian::Position::from_array(*pos.center_params(), rotation * *pos.as_array())
    }
}

/// Rotate a mean-equatorial J2000 position into `F` without requiring
/// `ReferenceCenter::Params = ()` (unlike `PositionAstroExt`).
fn rotate_mean_equatorial_position<C, U, Eph, F>(
    pos: &cartesian::Position<C, EquatorialMeanJ2000, U>,
    jd_reception: &JulianDate,
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
    let rotation = frame_rotation_with::<EquatorialMeanJ2000, F, _>(*jd_reception, ctx);
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

    /// Shared light-time observation geometry in sky frame `F`.
    fn observation_geometry<Body, E, F>(input: Self::Input) -> ObservationGeometry<Self, F>
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame;

    /// Body orientation in sky frame `F` with backend `E`.
    fn body_orientation<Body, E, F>(input: Self::Input) -> BodyOrientation
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame,
    {
        body_from_geometry(&Self::observation_geometry::<Body, E, F>(input))
    }

    /// Saturn ring orientation in sky frame `F` with backend `E`.
    fn saturn_ring_orientation<E, F>(input: Self::Input) -> SaturnRingOrientation
    where
        E: PlanetState<Saturn>,
        F: SkyOrientationFrame,
    {
        saturn_from_geometry(&Self::observation_geometry::<Saturn, E, F>(input))
    }
}

impl PlanetaryObservationCenter for Geocentric {
    type Input = JulianDate;

    fn observation_geometry<Body, E, F>(jd: Self::Input) -> ObservationGeometry<Self, F>
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let states = light_time_states::<E, Body>(&jd, &ctx, None);
        finish_geocentric_geometry::<Body, E, F>(&jd, &ctx, &states)
    }
}

impl PlanetaryObservationCenter for Topocentric {
    type Input = (JulianDate, Geodetic<ECEF>);

    fn observation_geometry<Body, E, F>((jd, site): Self::Input) -> ObservationGeometry<Self, F>
    where
        Body: HasIauRotation,
        E: PlanetState<Body>,
        F: SkyOrientationFrame,
    {
        let ctx = AstroContext::<E>::default();
        let site_ecl = site_geocentric_ecliptic_au(site, jd, &ctx);
        let states = light_time_states::<E, Body>(&jd, &ctx, Some(&site_ecl));
        finish_topocentric_geometry::<Body, E, F>(&jd, site, &ctx, &states)
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
    /// `None` when the planet–observer direction cannot be normalised.
    pub ring_opening_observer: Option<Radians>,
    /// Signed ring opening as seen from the Sun. Positive ⇒ north face
    /// illuminated. Independent of observation center and sky frame when
    /// defined. `None` when the heliocentric direction cannot be normalised.
    pub ring_opening_sun: Option<Radians>,
}

/// Shared sky-plane geometry for one observation (center `C`, frame `F`).
#[derive(Debug, Clone)]
pub struct ObservationGeometry<C: ReferenceCenter, F: ReferenceFrame> {
    /// Planet relative to the observation center.
    pub planet: cartesian::Position<C, F, AstronomicalUnit>,
    /// Sun relative to the observation center.
    pub sun: cartesian::Position<C, F, AstronomicalUnit>,
    /// IAU north pole in the sky frame (evaluated at the retarded epoch).
    pub pole: spherical::Direction<F>,
    /// Planet heliocentric vector (for solar opening / illumination geometry).
    pub planet_helio: cartesian::Position<Heliocentric, F, AstronomicalUnit>,
    /// Light-time retarded TT epoch used for the planet and IAU pole.
    pub emission_epoch: JulianDate,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Planetary disc orientation with compile-time body, backend, center, and frame.
///
/// Applies [`ORIENTATION_CORRECTIONS`] (light-time only). See the
/// [module documentation](self).
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
/// assert!(ori.ring_opening_observer.unwrap().to::<siderust::qtty::Degree>().abs() < Degrees::new(90.0));
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
/// # Examples
///
/// ```rust
/// use siderust::ephemeris::Vsop87Ephemeris;
/// use siderust::event::planetary::saturn_ring_orientation;
/// use siderust::qtty::Degrees;
///
/// let ori = saturn_ring_orientation::<Vsop87Ephemeris>(siderust::J2000);
/// assert!(ori.north_pole_position_angle.unwrap().is_finite());
/// assert!(ori.ring_opening_observer.unwrap().to::<siderust::qtty::Degree>().abs() < Degrees::new(90.0));
/// ```
#[inline]
pub fn saturn_ring_orientation<E: PlanetState<Saturn>>(jd: JulianDate) -> SaturnRingOrientation {
    saturn_rings::<E, Geocentric, EquatorialTrueOfDate>(jd)
}

impl Mars {
    /// Geocentric Mars orientation using [`DefaultEphemeris`] in
    /// [`EquatorialTrueOfDate`].
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
    planet_geo: GeoEcliptic,
    sun_geo: GeoEcliptic,
    planet_helio: HelioEcliptic,
    emission_epoch: JulianDate,
}

/// Geocentric site vector in mean ecliptic J2000 (AU), via the public
/// topocentric parallax path so EOP/ephemeris context stays consistent.
fn site_geocentric_ecliptic_au<E: Ephemeris>(
    site: Geodetic<ECEF>,
    jd: JulianDate,
    ctx: &AstroContext<E>,
) -> GeoEcliptic {
    let origin = GeoMean::new(0.0, 0.0, 0.0);
    let topo = to_topocentric_with(&origin, site, jd, ctx);
    // topo = origin − site = −site ⇒ site = −topo.
    let site_mean = GeoMean::new(-topo.x(), -topo.y(), -topo.z());
    TransformFrame::to_frame(&site_mean)
}

fn light_time_states<E, Body>(
    jd_reception: &JulianDate,
    _ctx: &AstroContext<E>,
    site_ecliptic: Option<&GeoEcliptic>,
) -> LightTimeStates
where
    E: PlanetState<Body>,
    Body: HasIauRotation,
{
    let earth = E::earth_heliocentric(*jd_reception);
    let observer = match site_ecliptic {
        None => HelioEcliptic::new(earth.x(), earth.y(), earth.z()),
        Some(site) => HelioEcliptic::new(
            earth.x() + site.x(),
            earth.y() + site.y(),
            earth.z() + site.z(),
        ),
    };

    let mut emission = *jd_reception;
    let mut planet = E::heliocentric(emission);
    let mut range = HelioEcliptic::new(
        planet.x() - observer.x(),
        planet.y() - observer.y(),
        planet.z() - observer.z(),
    );
    for _ in 0..LIGHT_TIME_ITERS {
        let tau_days = range.distance().value() / AU_PER_DAY_C.value();
        emission = JulianDate::new((jd_reception.raw() - Days::new(tau_days)).value());
        planet = E::heliocentric(emission);
        range = HelioEcliptic::new(
            planet.x() - observer.x(),
            planet.y() - observer.y(),
            planet.z() - observer.z(),
        );
    }

    let planet_geo = GeoEcliptic::new(
        planet.x() - earth.x(),
        planet.y() - earth.y(),
        planet.z() - earth.z(),
    );
    let sun_geo = GeoEcliptic::new(-earth.x(), -earth.y(), -earth.z());
    LightTimeStates {
        planet_geo,
        sun_geo,
        planet_helio: planet,
        emission_epoch: emission,
    }
}

fn finish_geocentric_geometry<Body, E, F>(
    jd_reception: &JulianDate,
    ctx: &AstroContext<E>,
    states: &LightTimeStates,
) -> ObservationGeometry<Geocentric, F>
where
    Body: HasIauRotation,
    E: Ephemeris,
    F: SkyOrientationFrame,
{
    let planet_mean: GeoMean = TransformFrame::to_frame(&states.planet_geo);
    let sun_mean: GeoMean = TransformFrame::to_frame(&states.sun_geo);
    let helio_mean: HelioMean = TransformFrame::to_frame(&states.planet_helio);
    ObservationGeometry {
        planet: F::position_from_mean_equatorial(&planet_mean, jd_reception, ctx),
        sun: F::position_from_mean_equatorial(&sun_mean, jd_reception, ctx),
        pole: iau_pole_in_frame::<Body, E, F>(&states.emission_epoch, jd_reception, ctx),
        planet_helio: F::position_from_mean_equatorial(&helio_mean, jd_reception, ctx),
        emission_epoch: states.emission_epoch,
    }
}

fn finish_topocentric_geometry<Body, E, F>(
    jd_reception: &JulianDate,
    site: Geodetic<ECEF>,
    ctx: &AstroContext<E>,
    states: &LightTimeStates,
) -> ObservationGeometry<Topocentric, F>
where
    Body: HasIauRotation,
    E: Ephemeris,
    F: SkyOrientationFrame,
{
    let planet_mean: GeoMean = TransformFrame::to_frame(&states.planet_geo);
    let sun_mean: GeoMean = TransformFrame::to_frame(&states.sun_geo);
    let helio_mean: HelioMean = TransformFrame::to_frame(&states.planet_helio);
    let planet_topo = to_topocentric_with(&planet_mean, site, *jd_reception, ctx);
    let sun_topo = to_topocentric_with(&sun_mean, site, *jd_reception, ctx);
    ObservationGeometry {
        planet: F::topocentric_from_mean_equatorial::<E, _>(&planet_topo, jd_reception, ctx),
        sun: F::topocentric_from_mean_equatorial::<E, _>(&sun_topo, jd_reception, ctx),
        pole: iau_pole_in_frame::<Body, E, F>(&states.emission_epoch, jd_reception, ctx),
        planet_helio: F::position_from_mean_equatorial(&helio_mean, jd_reception, ctx),
        emission_epoch: states.emission_epoch,
    }
}

fn iau_pole_in_frame<Body, E, F>(
    emission_epoch: &JulianDate,
    jd_reception: &JulianDate,
    ctx: &AstroContext<E>,
) -> spherical::Direction<F>
where
    Body: HasIauRotation,
    E: Ephemeris,
    F: SkyOrientationFrame,
{
    let params = &Body::ROTATION;
    let pole_icrs = spherical::Direction::<ICRS>::new(
        params.alpha0(*emission_epoch),
        params.delta0(*emission_epoch),
    );
    F::direction_from_icrs(&pole_icrs, jd_reception, ctx)
}

fn body_from_geometry<C, F>(geo: &ObservationGeometry<C, F>) -> BodyOrientation
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    BodyOrientation {
        north_pole_position_angle: bearing_to_pole(&geo.planet, &geo.pole),
        bright_limb_position_angle: bearing_to_position(&geo.planet, &geo.sun),
    }
}

fn saturn_from_geometry<C, F>(geo: &ObservationGeometry<C, F>) -> SaturnRingOrientation
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    SaturnRingOrientation {
        north_pole_position_angle: bearing_to_pole(&geo.planet, &geo.pole),
        ring_opening_observer: signed_opening(&geo.pole, &geo.planet),
        ring_opening_sun: signed_opening(&geo.pole, &geo.planet_helio),
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
) -> Option<Radians>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let dir = planet.direction()?;
    let pole_cart = pole.to_cartesian();
    // Planetocentric latitude of the observer: asin(n̂ · (−r̂)).
    let sin_b = -(pole_cart.dot(&dir));
    Some(Radians::new(clamp_unit(sin_b).asin()))
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
    fn orientation_policy_is_light_time_only() {
        assert_eq!(ORIENTATION_CORRECTIONS, CorrectionPolicy::LIGHT_TIME);
        assert!(!ORIENTATION_CORRECTIONS.contains(CorrectionPolicy::ABERRATION));
        assert!(!ORIENTATION_CORRECTIONS.contains(CorrectionPolicy::LIGHT_DEFLECTION));
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
        let b = ori.ring_opening_observer.unwrap().to::<Degree>();
        let bp = ori.ring_opening_sun.unwrap().to::<Degree>();
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
    fn geocentric_vs_topocentric_two_sites() {
        // Meeus 42.a epoch; independent check that parallax + observer light-time
        // shift the pole PA at the ~0.001°–0.1° level between distant sites.
        let jd = JulianDate::new(2_448_935.500_683);
        let geo = Mars::orientation(jd);
        let greenwich =
            Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
        let paranal =
            Geodetic::<ECEF>::new(Degrees::new(-70.4), Degrees::new(-24.6), Meters::new(0.0));
        let topo_g = Mars::orientation_at::<Topocentric, EquatorialTrueOfDate>((jd, greenwich));
        let topo_p = Mars::orientation_at::<Topocentric, EquatorialTrueOfDate>((jd, paranal));
        let g = geo.north_pole_position_angle.unwrap();
        let tg = topo_g.north_pole_position_angle.unwrap();
        let tp = topo_p.north_pole_position_angle.unwrap();
        assert!((g - tg).abs() > Radians::new(1e-8));
        assert!((g - tp).abs() > Radians::new(1e-8));
        assert!((tg - tp).abs() > Radians::new(1e-8));
        // Saturn openings: observer opening must react to site; solar opening must not.
        let sat_g =
            Saturn::ring_orientation_at::<Topocentric, EquatorialTrueOfDate>((jd, greenwich));
        let sat_p = Saturn::ring_orientation_at::<Topocentric, EquatorialTrueOfDate>((jd, paranal));
        assert!(
            (sat_g.ring_opening_observer.unwrap() - sat_p.ring_opening_observer.unwrap()).abs()
                > Radians::new(1e-10)
        );
        assert!(
            (sat_g.ring_opening_sun.unwrap() - sat_p.ring_opening_sun.unwrap()).abs()
                < Radians::new(1e-12)
        );
    }

    #[test]
    fn equatorial_frames_transform_and_openings_are_invariant() {
        let jd = JulianDate::new(2_448_972.500_68);
        let tod = saturn_rings::<DefaultEphemeris, Geocentric, EquatorialTrueOfDate>(jd);
        let mean = saturn_rings::<DefaultEphemeris, Geocentric, EquatorialMeanJ2000>(jd);
        let icrs = saturn_rings::<DefaultEphemeris, Geocentric, ICRS>(jd);
        assert!(
            (tod.ring_opening_observer.unwrap() - mean.ring_opening_observer.unwrap()).abs()
                < Radians::new(1e-10)
        );
        assert!(
            (tod.ring_opening_observer.unwrap() - icrs.ring_opening_observer.unwrap()).abs()
                < Radians::new(1e-10)
        );
        assert!(
            (tod.ring_opening_sun.unwrap() - mean.ring_opening_sun.unwrap()).abs()
                < Radians::new(1e-10)
        );
        assert!(
            (tod.ring_opening_sun.unwrap() - icrs.ring_opening_sun.unwrap()).abs()
                < Radians::new(1e-10)
        );
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
    fn retarded_pole_epoch_differs_from_reception_for_saturn() {
        let jd = JulianDate::new(2_448_972.500_68);
        let geo =
            Geocentric::observation_geometry::<Saturn, DefaultEphemeris, EquatorialTrueOfDate>(jd);
        assert!(
            (geo.emission_epoch.raw().value() - jd.raw().value()).abs() > 1e-6,
            "Saturn light-time should retard the emission epoch by minutes"
        );
        // Retardation for Saturn at ~10 AU is roughly 80 minutes ≈ 0.055 day.
        let dt_days = jd.raw().value() - geo.emission_epoch.raw().value();
        assert!(
            dt_days > 0.04 && dt_days < 0.08,
            "unexpected Saturn light-time {dt_days} d"
        );
    }

    #[test]
    fn signed_opening_none_for_zero_vector_and_zero_is_valid() {
        let pole = spherical::Direction::<EquatorialTrueOfDate>::new(
            Degrees::new(0.0),
            Degrees::new(90.0),
        );
        let zero = cartesian::Position::<Geocentric, EquatorialTrueOfDate, AstronomicalUnit>::new(
            0.0, 0.0, 0.0,
        );
        assert_eq!(signed_opening(&pole, &zero), None);

        // Planet along +X, pole along +Z ⇒ edge-on rings, opening = 0.
        let edge = cartesian::Position::<Geocentric, EquatorialTrueOfDate, AstronomicalUnit>::new(
            1.0, 0.0, 0.0,
        );
        let opening = signed_opening(&pole, &edge).unwrap();
        assert!(opening.abs() < Radians::new(1e-12));
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
            let b = ori.ring_opening_observer.unwrap();
            let bp = ori.ring_opening_sun.unwrap();
            assert!(b.abs() <= half_pi + Radians::new(1e-12));
            assert!(bp.abs() <= half_pi + Radians::new(1e-12));
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
            assert!(sat.ring_opening_observer.unwrap().is_finite());
            assert!(sat.ring_opening_sun.unwrap().is_finite());
        }
    }

    #[test]
    fn vsop87_backend_selection_matches_default_ephemeris_alias() {
        use crate::ephemeris::Vsop87Ephemeris;
        let jd = JulianDate::new(2_448_935.500_683);
        assert_eq!(
            mars_orientation::<Vsop87Ephemeris>(jd),
            mars_orientation::<DefaultEphemeris>(jd)
        );
        assert_eq!(
            saturn_ring_orientation::<Vsop87Ephemeris>(jd),
            saturn_ring_orientation::<DefaultEphemeris>(jd)
        );
    }
}
