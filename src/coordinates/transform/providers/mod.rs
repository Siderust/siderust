// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Transformation Provider Traits
//!
//! This module defines the provider traits for computing time-dependent
//! coordinate transformations. Providers are implemented for specific
//! frame/center pairs and return `affn` operators.
//!
//! ## Architecture
//!
//! The transformation system uses a hub-and-spoke model to avoid a
//! combinatorial explosion of implementations:
//!
//! - **Frame hub**: `ICRS` is the canonical inertial frame. Most frame
//!   rotations are composed through `ICRS`.
//! - **Center hub**: `Barycentric` is the canonical origin. Most center
//!   shifts are composed through `Barycentric`.
//!
//! ## Module Layout
//!
//! To keep this maintainable as new frames and centers are added, providers
//! are split by domain:
//!
//! - `frames_inertial`: ICRS, ICRF, J2000, EME2000, GCRS, mean/true equators.
//! - `frames_earth`: Earth orientation chain (`CIRS`, `TIRS`, `ITRF`, `ECEF`).
//! - `frames_catalog`: Galactic and FK4 catalog frames.
//! - `frames_teme`: TEME support for SGP4/TLE workflows.
//! - `frames_planetary`: IAU body-fixed planetary frames.
//! - `centers_standard`: barycentric, heliocentric and geocentric shifts.
//! - `centers_planetary`: planetocentric, plutocentric and selenocentric shifts.
//!
//! Shared composition helpers live in this root module so each submodule can
//! stay small and focused.

use crate::astro::earth_rotation::jd_ut1_from_tt_eop;
use crate::astro::earth_rotation_provider::nutation_with_celestial_pole_offsets;
use crate::astro::eop::EopProvider;
use crate::astro::nutation::NutationModel;
use crate::astro::{cio, era, polar_motion, precession, HasIauRotation, IauRotationParams};
use crate::calculus::ephemeris::Ephemeris;
use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::{
    Barycentric, Geocentric, Heliocentric, Jovicentric, Marscentric, Mercurycentric,
    Neptunocentric, Plutocentric, Saturnocentric, Selenocentric, Uranocentric, Venuscentric,
};
use crate::coordinates::frames::planetary::{
    JupiterSystemIII, MarsFixed, MercuryFixed, MoonPrincipalAxes, NeptuneFixed, PlutoFixed,
    SaturnFixed, UranusFixed, VenusFixed,
};
use crate::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate, Galactic,
    CIRS, ECEF, EME2000, FK4B1950, GCRS as GCRSFrame, ICRF, ICRS, ITRF, TEME, TIRS,
};
use crate::coordinates::transform::context::AstroContext;
use crate::time::JulianDate;
use affn::Rotation3;

mod centers_planetary;
mod centers_standard;
mod frames_catalog;
mod frames_earth;
mod frames_inertial;
mod frames_planetary;
mod frames_teme;

/// Trait for computing rotation matrices between reference frames.
///
/// Implementations provide the time-dependent rotation from frame `F1` to
/// frame `F2`. The rotation is computed at a given Julian Date using the
/// provided astronomical context.
///
/// # Type Parameters
///
/// - `F1`: Source reference frame.
/// - `F2`: Target reference frame.
///
/// # Sign Convention
///
/// The returned rotation transforms vectors from `F1` to `F2`:
/// ```text
/// v_F2 = rotation(F1 → F2) * v_F1
/// ```
pub trait FrameRotationProvider<F1, F2> {
    /// Computes the rotation matrix from frame `F1` to frame `F2`.
    fn rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3;
}

/// Trait for computing translation vectors between reference centers.
///
/// Implementations provide the time-dependent translation from center `C1`
/// to center `C2`, expressed in frame `F`. The result is a typed quantity
/// array in [`AstronomicalUnit`](qtty::AstronomicalUnit), matching the
/// canonical unit of the ephemeris providers.
pub trait CenterShiftProvider<C1, C2, F> {
    /// Computes the translation vector from center `C1` to center `C2`.
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [qtty::Quantity<qtty::AstronomicalUnit>; 3];
}

/// Identity rotation: same frame to same frame.
impl<F> FrameRotationProvider<F, F> for ()
where
    F: affn::ReferenceFrame,
{
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

/// Identity shift: same center to same center.
impl<C, F> CenterShiftProvider<C, C, F> for ()
where
    C: affn::ReferenceCenter,
    F: affn::ReferenceFrame,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [qtty::Quantity<qtty::AstronomicalUnit>; 3] {
        use qtty::AstronomicalUnit;
        [
            qtty::Quantity::<AstronomicalUnit>::new(0.0),
            qtty::Quantity::<AstronomicalUnit>::new(0.0),
            qtty::Quantity::<AstronomicalUnit>::new(0.0),
        ]
    }
}

/// Composes two frame rotations: `F1 → FM → F2`.
///
/// Equivalent to `R(FM→F2) · R(F1→FM)`, applied left-to-right.
#[inline]
fn compose_rotation<F1, FM, F2, Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Rotation3
where
    (): FrameRotationProvider<F1, FM>,
    (): FrameRotationProvider<FM, F2>,
{
    let left = <() as FrameRotationProvider<F1, FM>>::rotation(jd, ctx);
    let right = <() as FrameRotationProvider<FM, F2>>::rotation(jd, ctx);
    right * left
}

/// Returns `R(F2→F1)` by inverting `R(F1→F2)`.
#[inline]
fn inverse_rotation<F1, F2, Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Rotation3
where
    (): FrameRotationProvider<F2, F1>,
{
    <() as FrameRotationProvider<F2, F1>>::rotation(jd, ctx).inverse()
}

/// Convenience alias for a 3-component shift in astronomical units.
type AuShift = [qtty::Quantity<qtty::AstronomicalUnit>; 3];

/// Negates a 3-component shift vector (element-wise sign flip).
#[inline]
fn negate_shift(shift: AuShift) -> AuShift {
    [-shift[0], -shift[1], -shift[2]]
}

/// Adds two 3-component shift vectors element-wise.
#[inline]
fn add_shifts(lhs: AuShift, rhs: AuShift) -> AuShift {
    [lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]]
}

/// Returns `shift(C1→C2)` by negating `shift(C2→C1)`.
#[inline]
fn inverse_shift<C1, C2, F, Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> AuShift
where
    (): CenterShiftProvider<C2, C1, F>,
{
    negate_shift(<() as CenterShiftProvider<C2, C1, F>>::shift(jd, ctx))
}

/// Composes two center shifts: `C1 → CM → C2`.
#[inline]
fn compose_shift<C1, CM, C2, F, Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> AuShift
where
    (): CenterShiftProvider<C1, CM, F>,
    (): CenterShiftProvider<CM, C2, F>,
{
    let left = <() as CenterShiftProvider<C1, CM, F>>::shift(jd, ctx);
    let right = <() as CenterShiftProvider<CM, C2, F>>::shift(jd, ctx);
    add_shifts(left, right)
}

/// Direct GCRS → CIRS rotation from the IAU 2000/2006 CIO-based chain.
///
/// Applies IERS celestial-pole corrections (dX, dY) from the EOP provider
/// before computing the CIP/CIO-based rotation matrix.
#[inline]
fn gcrs_to_cirs_rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Rotation3 {
    let eop = ctx.eop_at(jd);
    let (dpsi, deps) = nutation_with_celestial_pole_offsets(jd, eop);
    let cip = cio::cip_cio(jd, dpsi, deps);
    cio::gcrs_to_cirs_matrix(cip.x, cip.y, cip.s)
}

/// Direct CIRS → TIRS rotation from the Earth Rotation Angle.
///
/// Converts the epoch from TT to UT1 using the EOP ΔUT1 correction,
/// then builds a single Rz rotation by the negative ERA.
#[inline]
fn cirs_to_tirs_rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Rotation3 {
    let eop = ctx.eop_at(jd);
    let jd_ut1 = jd_ut1_from_tt_eop(jd, &eop);
    Rotation3::rz(-era::earth_rotation_angle(jd_ut1))
}

/// Direct TIRS → ITRF rotation from polar motion.
///
/// Uses the IERS polar-motion parameters (xp, yp) from the EOP provider
/// to build the W(t) rotation matrix.
#[inline]
fn tirs_to_itrf_rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> Rotation3 {
    let eop = ctx.eop_at(jd);
    polar_motion::polar_motion_matrix_from_eop(eop.xp, eop.yp, jd)
}

/// Rotates a typed ecliptic position into target frame `F`, returning
/// a typed `[Quantity<AstronomicalUnit>; 3]` shift vector.
///
/// Uses [`Rotation3`]'s built-in `Mul<Position>` support to rotate the
/// position in typed quantity space — no raw `f64` extraction needed.
/// The result preserves the AU unit through the rotation, then extracts
/// the three components as typed quantities.
///
/// The center type `C` is irrelevant to the rotation and is only present
/// so callers can pass any ecliptic-plane position directly without
/// first converting it to a specific center.
#[inline]
fn rotate_shift_from_ecliptic<C, F, Eph, Eop: EopProvider, Nut: NutationModel>(
    pos: Position<C, EclipticMeanJ2000, qtty::AstronomicalUnit>,
    jd: JulianDate,
    ctx: &AstroContext<Eph, Eop, Nut>,
) -> AuShift
where
    C: affn::ReferenceCenter,
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    let rot = <() as FrameRotationProvider<EclipticMeanJ2000, F>>::rotation(jd, ctx);
    let rotated = rot * pos;
    [rotated.x(), rotated.y(), rotated.z()]
}

macro_rules! impl_via_icrs_bidirectional {
    ($src:ty, $dst:ty) => {
        impl FrameRotationProvider<$src, $dst> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                compose_rotation::<$src, ICRS, $dst, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl FrameRotationProvider<$dst, $src> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                inverse_rotation::<$dst, $src, Eph, Eop, Nut>(jd, ctx)
            }
        }
    };
}

pub(crate) use impl_via_icrs_bidirectional;

/// Implements the reverse and composed `CenterShiftProvider` impls for a
/// body center that already has a direct `Center → Barycentric` impl.
///
/// Given a primary `Center → Barycentric` implementation, this macro
/// generates the five remaining directional impls:
///
/// - `Barycentric → Center` (negate)
/// - `Heliocentric → Center` (compose via Barycentric)
/// - `Center → Heliocentric` (negate)
/// - `Geocentric → Center` (compose via Barycentric)
/// - `Center → Geocentric` (negate)
macro_rules! impl_reverse_center_shifts {
    ($center:ty) => {
        impl<F> CenterShiftProvider<Barycentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                inverse_shift::<Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<Heliocentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                compose_shift::<Heliocentric, Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<$center, Heliocentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                inverse_shift::<$center, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<Geocentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                compose_shift::<Geocentric, Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<$center, Geocentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                inverse_shift::<$center, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
            }
        }
    };
}

pub(crate) use impl_reverse_center_shifts;

/// Computes the body-fixed → ICRS rotation matrix from IAU rotation parameters.
///
/// Builds a ZXZ Euler rotation from the IAU pole direction (α₀, δ₀) and
/// prime-meridian angle (W). All angles are converted to radians via
/// `qtty` before being passed to [`Rotation3::from_euler_zxz`].
#[inline]
fn iau_body_fixed_to_icrs(params: &IauRotationParams, jd: JulianDate) -> Rotation3 {
    let alpha0 = params.alpha0(jd).to::<qtty::Radian>();
    let delta0 = params.delta0(jd).to::<qtty::Radian>();
    let w = params.w(jd).to::<qtty::Radian>();
    let right_angle = qtty::Degrees::new(90.0).to::<qtty::Radian>();

    let z1 = -w;
    let x = -(right_angle - delta0);
    let z2 = -(alpha0 + right_angle);
    Rotation3::from_euler_zxz(z1, x, z2)
}

/// Computes the center shift from `C1` to `C2` in frame `F`.
#[inline]
pub fn center_shift<C1, C2, F>(
    jd: JulianDate,
    ctx: &AstroContext,
) -> [qtty::Quantity<qtty::AstronomicalUnit>; 3]
where
    (): CenterShiftProvider<C1, C2, F>,
{
    <() as CenterShiftProvider<C1, C2, F>>::shift(jd, ctx)
}

/// Computes the rotation matrix from frame `F1` to frame `F2`.
#[inline]
pub fn frame_rotation<F1, F2>(jd: JulianDate, ctx: &AstroContext) -> Rotation3
where
    (): FrameRotationProvider<F1, F2>,
{
    <() as FrameRotationProvider<F1, F2>>::rotation(jd, ctx)
}

#[cfg(test)]
mod tests;
