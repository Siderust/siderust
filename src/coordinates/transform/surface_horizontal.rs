// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Generic body-surface geometric Alt/Az transforms
//!
//! ## Scientific scope
//!
//! Converts an inertial (or body-centred) target state into local horizontal
//! altitude and azimuth for an observer on a celestial body's reference
//! surface. The pipeline is body-agnostic:
//!
//! ```text
//! target (center C, frame F)
//!   → body-centred (B::Center) in ICRS
//!   → body-fixed (B::FixedFrame)
//!   → observer-relative line of sight
//!   → local ENU projection
//!   → Horizontal (alt, az)
//! ```
//!
//! Body orientation uses the body's IAU WGCCRE rotation model at the
//! requested TT epoch. Earth sidereal-time (GMST/GAST) logic is **not**
//! reused for the Moon, Mars, or other bodies.
//!
//! Results are **geometric**: no refraction, aberration, light deflection,
//! light-time, or terrain masking.
//!
//! ## Technical scope
//!
//! - [`position_to_horizontal`] — finite-distance target with topocentric parallax
//! - [`direction_to_horizontal`] — direction-only (no invented parallax)
//! - [`horizontal_to_direction`] — inverse: local horizontal → body-fixed /
//!   inertial direction
//! - [`SurfaceHorizontalExt`] — ergonomic method on [`SurfaceSite`]
//!
//! Azimuth convention: north = 0°, east = 90°, increasing north through east.
//! Altitude: geometric horizon = 0°, zenith = +90°.
//!
//! ## References
//!
//! - IAU WGCCRE (2015/2018) body-fixed frames
//! - Vallado (2013) §3.4 (ENU / az-el geometry)
//! - See [`crate::coordinates::transform::horizontal`] for the Earth-only
//!   GAST equatorial path retained for backward compatibility

use crate::coordinates::cartesian::{Direction, Position};
use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::{Horizontal, ReferenceFrame, ICRS};
use crate::coordinates::spherical;
use crate::coordinates::surface::{
    BodySurfaceModel, SurfaceError, SurfaceObserver, SurfaceObserverParams, SurfaceSite,
};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::TransformContext;
use crate::coordinates::transform::ext::{DirectionAstroExt, PositionAstroExt};
use crate::coordinates::transform::providers::{CenterShiftProvider, FrameRotationProvider};
use crate::ephemeris::Ephemeris;
use crate::qtty::{Degrees, LengthUnit, Meter, Quantity};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};
use crate::time::JulianDate;
use core::f64::consts::{PI, TAU};

/// Geometric horizontal observation from a surface site.
///
/// Carries the observer site in [`SurfaceObserver`] centre parameters so the
/// result cannot be silently reinterpreted at a different site.
#[derive(Debug, Clone, Copy)]
pub struct HorizontalObservation<B: BodySurfaceModel, U: LengthUnit = Meter> {
    /// Topocentric horizontal position (range preserved for finite targets).
    pub position: Position<SurfaceObserver<B>, Horizontal, U>,
}

impl<B: BodySurfaceModel, U: LengthUnit> HorizontalObservation<B, U> {
    /// Altitude above the geometric horizon.
    pub fn alt(&self) -> Degrees {
        let sph = affn::Position::<SurfaceObserver<B>, Horizontal, U>::to_spherical(&self.position);
        sph.alt()
    }

    /// Azimuth from north through east, normalised to `[0°, 360°)`.
    pub fn az(&self) -> Degrees {
        let sph = affn::Position::<SurfaceObserver<B>, Horizontal, U>::to_spherical(&self.position);
        sph.az()
    }

    /// Slant range from observer to target.
    pub fn range(&self) -> Quantity<U> {
        self.position.distance()
    }

    /// Observer site parameters attached to this observation.
    pub fn site_params(&self) -> &SurfaceObserverParams<B> {
        self.position.center_params()
    }

    /// Spherical horizontal direction (altitude, azimuth).
    pub fn direction_spherical(&self) -> spherical::Direction<Horizontal> {
        let sph = affn::Position::<SurfaceObserver<B>, Horizontal, U>::to_spherical(&self.position);
        sph.direction()
    }

    /// Cartesian horizontal direction.
    pub fn direction(&self) -> Direction<Horizontal> {
        self.direction_spherical().to_cartesian()
    }
}

/// Project a unit line-of-sight in a body-fixed frame onto local ENU / Horizontal.
fn los_to_horizontal(
    los_unit: [f64; 3],
    east: [f64; 3],
    north: [f64; 3],
    up: [f64; 3],
) -> Result<(Degrees, Degrees), SurfaceError> {
    let e = dot(los_unit, east);
    let n = dot(los_unit, north);
    let u = dot(los_unit, up);
    if ![e, n, u].iter().all(|v| v.is_finite()) {
        return Err(SurfaceError::NonFiniteInput {
            field: "line_of_sight",
        });
    }
    let horiz = (e * e + n * n).sqrt();
    let az = if horiz == 0.0 {
        // Zenith / nadir: azimuth undefined; adopt 0° by convention.
        0.0
    } else {
        e.atan2(n).rem_euclid(TAU)
    };
    let alt = u.atan2(horiz);
    Ok((Degrees::new(rad_to_deg(alt)), Degrees::new(rad_to_deg(az))))
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

/// Finite-distance geometric Alt/Az from a surface observer.
///
/// The target is translated to the body centre and rotated into the body's
/// fixed frame at `jd` (TT). The observer-relative vector includes
/// topocentric parallax. Direction-only inputs must use
/// [`direction_to_horizontal`] instead — this function will not invent a
/// distance.
///
/// # Type constraints
///
/// - `C` must admit a centre shift to `B::Center` in frame `F`
/// - `F` must rotate to [`ICRS`] (hub), and `ICRS` must rotate to `B::FixedFrame`
pub fn position_to_horizontal<B, C, F, U, Ctx>(
    target: &Position<C, F, U>,
    site: &SurfaceSite<B>,
    jd: &JulianDate,
    ctx: &Ctx,
) -> Result<HorizontalObservation<B, U>, SurfaceError>
where
    B: BodySurfaceModel,
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: From<Quantity<Meter>> + Into<Quantity<Meter>>,
    Ctx: TransformContext,
    Ctx::Eph: Ephemeris,
    (): CenterShiftProvider<C, B::Center, F>,
    (): FrameRotationProvider<F, ICRS>,
    (): FrameRotationProvider<ICRS, B::FixedFrame>,
    (): FrameRotationProvider<B::FixedFrame, ICRS>,
{
    B::ensure_supported()?;
    let params = site.to_params();

    // 1. Centre → body centre in source frame, then rotate to ICRS, then to body-fixed.
    let body_centred: Position<B::Center, F, U> =
        TransformCenter::to_center_with(target, (), *jd, ctx);
    let body_icrs: Position<B::Center, ICRS, U> = body_centred.to_frame_with(jd, ctx);
    let body_fixed: Position<B::Center, B::FixedFrame, U> = body_icrs.to_frame_with(jd, ctx);

    // 2. Observer in body-fixed metres, converted to U.
    let obs = site.observer_position::<U>()?;
    let dx = (body_fixed.x() - obs.x()).value();
    let dy = (body_fixed.y() - obs.y()).value();
    let dz = (body_fixed.z() - obs.z()).value();
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    if !r.is_finite() {
        return Err(SurfaceError::NonFiniteInput { field: "range" });
    }
    // Compare against a physical distance, never a raw numeric value in U:
    // 1e-6 AU would otherwise reject perfectly valid nearby targets.
    // Sub-micron residual after frame round-trips counts as coincident.
    let range_m: Quantity<Meter> = Quantity::<U>::new(r).into();
    if range_m.value() < 1e-6 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "zero observer–target range",
        });
    }
    let los = [dx / r, dy / r, dz / r];

    let enu = site.enu_basis()?;
    let (alt, az) = los_to_horizontal(los, enu.east, enu.north, enu.up)?;

    // Horizontal Cartesian: x=North, y=West, z=Zenith (affn Horizontal convention
    // used by existing Topocentric transforms).
    let az_r = deg_to_rad(az);
    let alt_r = deg_to_rad(alt);
    let hx = Quantity::<U>::new(r * alt_r.cos() * az_r.cos());
    let hy = Quantity::<U>::new(-(r * alt_r.cos() * az_r.sin()));
    let hz = Quantity::<U>::new(r * alt_r.sin());

    let position = Position::<SurfaceObserver<B>, Horizontal, U>::from_array(params, [hx, hy, hz]);
    Ok(HorizontalObservation { position })
}

/// Direction-only geometric Alt/Az (no topocentric parallax).
///
/// The input direction is rotated into the body-fixed frame and projected
/// onto the local ENU basis. No finite distance is assumed.
pub fn direction_to_horizontal<B, F, Ctx>(
    target: &Direction<F>,
    site: &SurfaceSite<B>,
    jd: &JulianDate,
    ctx: &Ctx,
) -> Result<Direction<Horizontal>, SurfaceError>
where
    B: BodySurfaceModel,
    F: ReferenceFrame,
    Ctx: TransformContext,
    (): FrameRotationProvider<F, ICRS>,
    (): FrameRotationProvider<ICRS, B::FixedFrame>,
{
    B::ensure_supported()?;
    let icrs: Direction<ICRS> = target.to_frame_with(jd, ctx);
    let fixed: Direction<B::FixedFrame> = icrs.to_frame_with(jd, ctx);
    let los = [fixed.x(), fixed.y(), fixed.z()];
    let n = norm(los);
    if !n.is_finite() || n == 0.0 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "zero-length direction",
        });
    }
    let los = [los[0] / n, los[1] / n, los[2] / n];
    let enu = site.enu_basis()?;
    let (alt, az) = los_to_horizontal(los, enu.east, enu.north, enu.up)?;
    let sph = spherical::Direction::<Horizontal>::new(alt, az);
    Ok(sph.to_cartesian())
}

/// Inverse: horizontal direction → body-fixed unit direction.
pub fn horizontal_to_body_fixed_direction<B>(
    horizontal: &Direction<Horizontal>,
    site: &SurfaceSite<B>,
) -> Result<Direction<B::FixedFrame>, SurfaceError>
where
    B: BodySurfaceModel,
{
    B::ensure_supported()?;
    let sph = horizontal.to_spherical();
    let alt = deg_to_rad(sph.alt());
    let az = deg_to_rad(sph.az());
    // Horizontal cartesian: x=N, y=W, z=Up → ENU components
    // e = sin(az) cos(alt), n = cos(az) cos(alt), u = sin(alt)
    let e = az.sin() * alt.cos();
    let n = az.cos() * alt.cos();
    let u = alt.sin();
    let enu = site.enu_basis()?;
    let mut v = [
        e * enu.east[0] + n * enu.north[0] + u * enu.up[0],
        e * enu.east[1] + n * enu.north[1] + u * enu.up[1],
        e * enu.east[2] + n * enu.north[2] + u * enu.up[2],
    ];
    let nv = norm(v);
    if !nv.is_finite() || nv == 0.0 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "zero inverse direction",
        });
    }
    v = [v[0] / nv, v[1] / nv, v[2] / nv];
    Ok(Direction::<B::FixedFrame>::new_unchecked(v[0], v[1], v[2]))
}

/// Inverse: horizontal direction → inertial/`F` direction at epoch `jd`.
pub fn horizontal_to_direction<B, F, Ctx>(
    horizontal: &Direction<Horizontal>,
    site: &SurfaceSite<B>,
    jd: &JulianDate,
    ctx: &Ctx,
) -> Result<Direction<F>, SurfaceError>
where
    B: BodySurfaceModel,
    F: ReferenceFrame,
    Ctx: TransformContext,
    (): FrameRotationProvider<B::FixedFrame, ICRS>,
    (): FrameRotationProvider<ICRS, F>,
{
    let fixed = horizontal_to_body_fixed_direction(horizontal, site)?;
    let icrs: Direction<ICRS> = fixed.to_frame_with(jd, ctx);
    Ok(icrs.to_frame_with(jd, ctx))
}

/// Ergonomic surface-site methods for geometric Alt/Az.
pub trait SurfaceHorizontalExt<B: BodySurfaceModel> {
    /// Finite-distance geometric Alt/Az for `target` at TT epoch `jd`.
    fn observe_position<C, F, U, Ctx>(
        &self,
        target: &Position<C, F, U>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> Result<HorizontalObservation<B, U>, SurfaceError>
    where
        C: ReferenceCenter<Params = ()>,
        F: ReferenceFrame,
        U: LengthUnit,
        Quantity<U>: From<Quantity<Meter>> + Into<Quantity<Meter>>,
        Ctx: TransformContext,
        Ctx::Eph: Ephemeris,
        (): CenterShiftProvider<C, B::Center, F>,
        (): FrameRotationProvider<F, ICRS>,
        (): FrameRotationProvider<ICRS, B::FixedFrame>,
        (): FrameRotationProvider<B::FixedFrame, ICRS>;

    /// Direction-only geometric Alt/Az (no parallax).
    fn observe_direction<F, Ctx>(
        &self,
        target: &Direction<F>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> Result<Direction<Horizontal>, SurfaceError>
    where
        F: ReferenceFrame,
        Ctx: TransformContext,
        (): FrameRotationProvider<F, ICRS>,
        (): FrameRotationProvider<ICRS, B::FixedFrame>;
}

impl<B: BodySurfaceModel> SurfaceHorizontalExt<B> for SurfaceSite<B> {
    fn observe_position<C, F, U, Ctx>(
        &self,
        target: &Position<C, F, U>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> Result<HorizontalObservation<B, U>, SurfaceError>
    where
        C: ReferenceCenter<Params = ()>,
        F: ReferenceFrame,
        U: LengthUnit,
        Quantity<U>: From<Quantity<Meter>> + Into<Quantity<Meter>>,
        Ctx: TransformContext,
        Ctx::Eph: Ephemeris,
        (): CenterShiftProvider<C, B::Center, F>,
        (): FrameRotationProvider<F, ICRS>,
        (): FrameRotationProvider<ICRS, B::FixedFrame>,
        (): FrameRotationProvider<B::FixedFrame, ICRS>,
    {
        position_to_horizontal(target, self, jd, ctx)
    }

    fn observe_direction<F, Ctx>(
        &self,
        target: &Direction<F>,
        jd: &JulianDate,
        ctx: &Ctx,
    ) -> Result<Direction<Horizontal>, SurfaceError>
    where
        F: ReferenceFrame,
        Ctx: TransformContext,
        (): FrameRotationProvider<F, ICRS>,
        (): FrameRotationProvider<ICRS, B::FixedFrame>,
    {
        direction_to_horizontal(target, self, jd, ctx)
    }
}

#[inline]
fn deg_to_rad(d: Degrees) -> f64 {
    d.value() * (PI / 180.0)
}

#[inline]
fn rad_to_deg(r: f64) -> f64 {
    r * (180.0 / PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::cartesian::Position;
    use crate::coordinates::centers::Selenocentric;
    use crate::coordinates::frames::MoonPrincipalAxes;
    use crate::coordinates::surface::MoonSurface;
    use crate::coordinates::transform::AstroContext;
    use crate::qtty::{Kilometer, DEG, KM, M};

    #[test]
    fn lunar_zenith_is_altitude_90() {
        let site = SurfaceSite::<MoonSurface>::new(10.0 * DEG, 20.0 * DEG, 0.0 * M).unwrap();
        let [x, y, z] = site.body_fixed_cartesian().unwrap();
        // Target further along the same radius vector → zenith.
        let scale = 2.0;
        let target = Position::<Selenocentric, MoonPrincipalAxes, Kilometer>::new(
            Quantity::<Kilometer>::new(x.value() * scale / 1000.0),
            Quantity::<Kilometer>::new(y.value() * scale / 1000.0),
            Quantity::<Kilometer>::new(z.value() * scale / 1000.0),
        );
        let jd = crate::J2000;
        let ctx =
            AstroContext::<crate::coordinates::transform::context::DefaultEphemeris>::default();
        let obs = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
        assert!((obs.alt().value() - 90.0).abs() < 1e-6, "alt={}", obs.alt());
        let expected_range_km = 1_737.4;
        assert!(
            (obs.range().value() - expected_range_km).abs() < 1e-3,
            "range={}",
            obs.range()
        );
        let _ = KM;
    }

    #[test]
    fn direction_roundtrip_lunar() {
        let site = SurfaceSite::<MoonSurface>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        let jd = crate::J2000;
        let ctx =
            AstroContext::<crate::coordinates::transform::context::DefaultEphemeris>::default();
        // Local east in body-fixed at (0,0): east = (0,1,0)
        let east_fixed = Direction::<MoonPrincipalAxes>::new_unchecked(0.0, 1.0, 0.0);
        let icrs: Direction<ICRS> = east_fixed.to_frame_with(&jd, &ctx);
        let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
        let sph = horiz.to_spherical();
        assert!((sph.alt().value()).abs() < 1e-6, "alt={}", sph.alt());
        let az = sph.az().value();
        assert!((az - 90.0).abs() < 1e-4, "az={az}");

        let back: Direction<ICRS> = horizontal_to_direction(&horiz, &site, &jd, &ctx).unwrap();
        let ang = icrs.angle_to(&back).value();
        assert!(ang < 1e-10, "roundtrip angle={ang}");
    }
}
