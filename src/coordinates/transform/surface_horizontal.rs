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
//! Altitude: geometric horizon = 0°, zenith = +90°. At zenith/nadir azimuth is
//! physically undefined and is reported as `0°` by convention whenever the
//! horizontal ENU projection falls below a scale-free singularity floor tied
//! to `f64` round-off and the ENU orthonormality tolerance (`1e-12`).
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
    /// Observer-tagged topocentric horizontal position.
    position: Position<SurfaceObserver<B>, Horizontal, U>,
}

impl<B: BodySurfaceModel, U: LengthUnit> HorizontalObservation<B, U> {
    // Runtime frame callers must validate their frame/site association first.
    #[cfg(feature = "spice")]
    pub(crate) fn from_fixed_relative(
        site: &SurfaceSite<B>,
        relative: [Quantity<U>; 3],
    ) -> Result<Self, SurfaceError> {
        let v = relative.map(|q| q.value());
        if v.iter().any(|x| !x.is_finite()) {
            return Err(SurfaceError::NonFiniteInput {
                field: "relative position",
            });
        }
        if !norm(v).is_finite() {
            return Err(SurfaceError::NonFiniteInput { field: "range" });
        }
        if norm(v) <= 0.0 {
            return Err(SurfaceError::DegenerateGeometry {
                reason: "coincident target and site",
            });
        }
        let basis = site.enu_basis()?;
        let horizontal = [dot(v, basis.north), dot(v, basis.east), dot(v, basis.up)];
        if horizontal.iter().any(|v| !v.is_finite()) {
            return Err(SurfaceError::NonFiniteInput {
                field: "horizontal position",
            });
        }
        Ok(Self {
            position: Position::from_array(site.to_params(), horizontal.map(Quantity::new)),
        })
    }

    /// The validated observer-tagged Cartesian position.
    pub fn position(&self) -> &Position<SurfaceObserver<B>, Horizontal, U> {
        &self.position
    }

    // affn's generic spherical conversion treats radii smaller than
    // f64::EPSILON in the *input unit* as zero. At micrometre-scale ranges
    // represented in AU that would silently erase a valid sky direction.
    // Use a scale-free Cartesian angular extraction instead.
    fn stable_angles(&self) -> (Degrees, Degrees) {
        let x = self.position.x().value();
        let y = self.position.y().value();
        let z = self.position.z().value();
        angles_from_horizontal_cartesian(x, y, z)
    }

    /// Altitude above the geometric horizon.
    pub fn alt(&self) -> Degrees {
        self.stable_angles().0
    }

    /// Azimuth from north through east, normalised to `[0°, 360°)`.
    pub fn az(&self) -> Degrees {
        self.stable_angles().1
    }

    /// Slant range from observer to target.
    pub fn range(&self) -> Quantity<U> {
        Quantity::new(norm([
            self.position.x().value(),
            self.position.y().value(),
            self.position.z().value(),
        ]))
    }

    /// Observer site parameters attached to this observation.
    pub fn site_params(&self) -> &SurfaceObserverParams<B> {
        self.position.center_params()
    }

    /// Spherical horizontal direction (altitude, azimuth).
    pub fn direction_spherical(&self) -> spherical::Direction<Horizontal> {
        let (altitude, azimuth) = self.stable_angles();
        spherical::Direction::<Horizontal>::new_unchecked(altitude, azimuth)
    }

    /// Cartesian horizontal direction.
    pub fn direction(&self) -> Direction<Horizontal> {
        self.direction_spherical().to_cartesian()
    }
}

/// Local geometric sky direction tied to its specific body and observer site.
///
/// The raw `Direction<Horizontal>` alone has no centre parameters; this
/// wrapper preserves them and prevents silently using the direction at another
/// body or location.
#[derive(Debug, Clone, Copy)]
pub struct HorizontalDirectionObservation<B: BodySurfaceModel> {
    /// Underlying horizontal unit direction.
    pub direction: Direction<Horizontal>,
    site_params: SurfaceObserverParams<B>,
}

impl<B: BodySurfaceModel> HorizontalDirectionObservation<B> {
    /// Observer parameters used to calculate this direction.
    pub fn site_params(&self) -> &SurfaceObserverParams<B> {
        &self.site_params
    }

    /// Horizontal sky direction in spherical coordinates.
    ///
    /// Uses the same singularity-aware Alt/Az extraction as
    /// [`HorizontalObservation`] so zenith/nadir azimuth stays at the
    /// documented `0°` convention.
    pub fn to_spherical(&self) -> spherical::Direction<Horizontal> {
        let (alt, az) = angles_from_horizontal_cartesian(
            self.direction.x(),
            self.direction.y(),
            self.direction.z(),
        );
        spherical::Direction::<Horizontal>::new_unchecked(alt, az)
    }

    /// Geometric altitude at the observer site.
    pub fn alt(&self) -> Degrees {
        angles_from_horizontal_cartesian(self.direction.x(), self.direction.y(), self.direction.z())
            .0
    }

    /// Azimuth (north through east) at the observer site.
    pub fn az(&self) -> Degrees {
        angles_from_horizontal_cartesian(self.direction.x(), self.direction.y(), self.direction.z())
            .1
    }
}

/// Inverse of a site-bound direction, rejecting attempts to reuse it at a
/// different observer site; the body itself is enforced by the type `B`.
pub fn horizontal_observation_to_direction<B, F, Ctx>(
    observation: &HorizontalDirectionObservation<B>,
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
    if observation.site_params != site.to_params() {
        return Err(SurfaceError::MismatchedObserver);
    }
    horizontal_to_direction(&observation.direction, site, jd, ctx)
}

/// Return a site- and body-bound horizontal sky direction instead of the
/// untyped `Direction<Horizontal>` returned by the low-level conversion.
pub fn observed_direction_to_horizontal<B, F, Ctx>(
    target: &Direction<F>,
    site: &SurfaceSite<B>,
    jd: &JulianDate,
    ctx: &Ctx,
) -> Result<HorizontalDirectionObservation<B>, SurfaceError>
where
    B: BodySurfaceModel,
    F: ReferenceFrame,
    Ctx: TransformContext,
    (): FrameRotationProvider<F, ICRS>,
    (): FrameRotationProvider<ICRS, B::FixedFrame>,
{
    let direction = direction_to_horizontal(target, site, jd, ctx)?;
    Ok(HorizontalDirectionObservation {
        direction,
        site_params: site.to_params(),
    })
}

/// Relative floor below which a horizontal ENU projection is treated as
/// singular for azimuth.
///
/// Two contributions, both scale-free via multiplication by `hypot(e,n,u)`:
///
/// 1. `32 ε` — local fused multiply–add noise when forming the ENU dots.
/// 2. `1e-12` — the orthonormality tolerance enforced by
///    [`crate::coordinates::surface::EnuBasis::assert_orthonormal`]. A
///    horizontal component below that level is not meaningful relative to
///    the basis itself, and also covers the ~1e-14 unit-LOS residual seen
///    when a short vertical offset is formed via ICRS subtraction of two
///    body-centred positions of order the body's radius.
///
/// This does **not** claim mathematical exactness from arbitrary floating
/// point. Near-zenith directions with `horiz ≳ 1e-10` (≳ 20 µas) still
/// retain a well-defined azimuth. Coarser “practically vertical” policies
/// remain the caller's responsibility.
const AZIMUTH_SINGULARITY_REL: f64 = {
    let eps_floor = 32.0 * f64::EPSILON;
    if eps_floor > 1e-12 {
        eps_floor
    } else {
        1e-12
    }
};

/// Azimuth from north-through-east ENU components, with a scale-free
/// zenith/nadir singularity.
///
/// Returns the documented convention `0` when
/// `sqrt(e² + n²) ≤ `[`AZIMUTH_SINGULARITY_REL`]` · hypot(e, n, u)`.
#[inline]
fn azimuth_from_enu(east: f64, north: f64, up: f64) -> f64 {
    let horiz = (east * east + north * north).sqrt();
    let radial = (east * east + north * north + up * up).sqrt();
    if !horiz.is_finite() || !radial.is_finite() || radial == 0.0 {
        return 0.0;
    }
    if horiz <= AZIMUTH_SINGULARITY_REL * radial {
        0.0
    } else {
        east.atan2(north).rem_euclid(TAU)
    }
}

/// Scale-free Alt/Az from affn Horizontal Cartesian `(N, E, Up)`.
#[inline]
fn angles_from_horizontal_cartesian(north: f64, east: f64, up: f64) -> (Degrees, Degrees) {
    let scale = north.abs().max(east.abs()).max(up.abs());
    if scale == 0.0 || !scale.is_finite() {
        return (Degrees::new(0.0), Degrees::new(0.0));
    }
    let (n, e, u) = (north / scale, east / scale, up / scale);
    let horiz = (n * n + e * e).sqrt();
    let altitude = rad_to_deg(u.atan2(horiz));
    let azimuth = rad_to_deg(azimuth_from_enu(e, n, u));
    (Degrees::new(altitude), Degrees::new(azimuth))
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
    let az = azimuth_from_enu(e, n, u);
    let alt = u.atan2(horiz);
    Ok((Degrees::new(rad_to_deg(alt)), Degrees::new(rad_to_deg(az))))
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    let scale = v[0].abs().max(v[1].abs()).max(v[2].abs());
    if scale == 0.0 {
        return 0.0;
    }
    let scaled = [v[0] / scale, v[1] / scale, v[2] / scale];
    scale * dot(scaled, scaled).sqrt()
}

/// Finite-distance geometric Alt/Az from a surface observer.
///
/// The target is translated to the body centre and rotated into the body's
/// fixed frame at `jd` (TT). The observer-relative vector includes
/// topocentric parallax. Direction-only inputs must use
/// [`direction_to_horizontal`] instead — this function will not invent a
/// distance. Coincident or numerically unresolved sightlines, and ranges
/// not representable in the requested output unit, return an error.
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
    // SurfaceSite is validated at construction; custom reference surfaces
    // are permitted for body models without a default surface.
    let params = site.to_params();

    // Canonicalise to metres BEFORE rotations, centre shifts and subtracting
    // observer coordinates. Working in arbitrary U (especially AU) makes a
    // near-surface range depend on the unit's floating-point spacing.
    let target_m: Position<C, F, Meter> = target.to_unit::<Meter>();
    let body_centred: Position<B::Center, F, Meter> =
        TransformCenter::to_center_with(&target_m, (), *jd, ctx);
    let body_icrs: Position<B::Center, ICRS, Meter> = body_centred.to_frame_with(jd, ctx);
    // Form the observer-relative vector in ICRS, then rotate only the unit
    // LOS into the body-fixed frame. Subtracting after a Fixed↔ICRS round-trip
    // of the absolute positions would leave an O(ε·R) residual that swamps a
    // short vertical offset and invents a spurious azimuth at zenith/nadir.
    let observer = site.observer_position::<Meter>()?;
    let observer_icrs: Position<B::Center, ICRS, Meter> = observer.to_frame_with(jd, ctx);

    let dx = (body_icrs.x() - observer_icrs.x()).value();
    let dy = (body_icrs.y() - observer_icrs.y()).value();
    let dz = (body_icrs.z() - observer_icrs.z()).value();
    if ![dx, dy, dz].iter().all(|value| value.is_finite()) {
        return Err(SurfaceError::NonFiniteInput {
            field: "observer_relative_position",
        });
    }
    // Scale before evaluating the norm: no overflow in squaring huge values
    // and no underflow when evaluating a tiny displacement.
    let scale = dx.abs().max(dy.abs()).max(dz.abs());
    let r_m = if scale == 0.0 {
        0.0
    } else {
        let scaled = [dx / scale, dy / scale, dz / scale];
        scale * (scaled[0] * scaled[0] + scaled[1] * scaled[1] + scaled[2] * scaled[2]).sqrt()
    };
    if !r_m.is_finite() {
        return Err(SurfaceError::NonFiniteInput { field: "range" });
    }
    // The tolerance is physical (metres), independent of requested U.
    // It also includes the round-off floor for subtracting two absolute
    // positions of order the body's radius.
    let absolute_scale = body_icrs
        .x()
        .value()
        .abs()
        .max(body_icrs.y().value().abs())
        .max(body_icrs.z().value().abs())
        .max(observer_icrs.x().value().abs())
        .max(observer_icrs.y().value().abs())
        .max(observer_icrs.z().value().abs());
    let resolution_m = (64.0 * f64::EPSILON * absolute_scale).max(1.0e-6);
    if r_m <= resolution_m {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "coincident or numerically indistinguishable observer and target",
        });
    }
    let los_icrs = Direction::<ICRS>::new_unchecked(dx / r_m, dy / r_m, dz / r_m);
    let los_fixed: Direction<B::FixedFrame> = los_icrs.to_frame_with(jd, ctx);
    let los = [los_fixed.x(), los_fixed.y(), los_fixed.z()];

    let enu = site.enu_basis()?;
    let (alt, az) = los_to_horizontal(los, enu.east, enu.north, enu.up)?;

    // Horizontal Cartesian: x=North, y=East, z=Zenith (affn spherical convention).
    let az_r = deg_to_rad(az);
    let alt_r = deg_to_rad(alt);
    let range = Quantity::<U>::from(Quantity::<Meter>::new(r_m));
    if !range.value().is_finite() {
        return Err(SurfaceError::NonFiniteInput {
            field: "range in output unit",
        });
    }
    if range.value() == 0.0 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "range is not representable in output unit",
        });
    }
    let hx = range * alt_r.cos() * az_r.cos();
    let hy = range * alt_r.cos() * az_r.sin();
    let hz = range * alt_r.sin();

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
    // Read canonical Horizontal axes directly: N, E, Up. Spherical
    // conversion would assign (alt, az)=(0, 0) to an unchecked zero vector.
    let length = norm([horizontal.x(), horizontal.y(), horizontal.z()]);
    if !length.is_finite() || length == 0.0 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "zero or non-finite horizontal direction",
        });
    }
    let n = horizontal.x() / length;
    let e = horizontal.y() / length;
    let u = horizontal.z() / length;
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
    ) -> Result<HorizontalDirectionObservation<B>, SurfaceError>
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
    ) -> Result<HorizontalDirectionObservation<B>, SurfaceError>
    where
        F: ReferenceFrame,
        Ctx: TransformContext,
        (): FrameRotationProvider<F, ICRS>,
        (): FrameRotationProvider<ICRS, B::FixedFrame>,
    {
        observed_direction_to_horizontal(target, self, jd, ctx)
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
    use crate::bodies::solar_system::Moon;
    use crate::coordinates::cartesian::Position;
    use crate::coordinates::centers::Selenocentric;
    use crate::coordinates::frames::MoonIauMeanEarth;
    use crate::coordinates::transform::AstroContext;
    use crate::qtty::{Kilometer, DEG, KM, M};

    #[test]
    fn lunar_zenith_is_altitude_90() {
        let site = SurfaceSite::<Moon>::new(10.0 * DEG, 20.0 * DEG, 0.0 * M).unwrap();
        let [x, y, z] = site.body_fixed_cartesian().unwrap();
        // Target further along the same radius vector → zenith.
        let scale = 2.0;
        let target = Position::<Selenocentric, MoonIauMeanEarth, Kilometer>::new(
            Quantity::<Kilometer>::new(x.value() * scale / 1000.0),
            Quantity::<Kilometer>::new(y.value() * scale / 1000.0),
            Quantity::<Kilometer>::new(z.value() * scale / 1000.0),
        );
        let jd = crate::J2000;
        let ctx =
            AstroContext::<crate::coordinates::transform::context::DefaultEphemeris>::default();
        let obs = position_to_horizontal(&target, &site, &jd, &ctx).unwrap();
        assert!((obs.alt().value() - 90.0).abs() < 1e-6, "alt={}", obs.alt());
        assert_eq!(obs.az().value(), 0.0, "zenith az convention");
        let expected_range_km = 1_737.4;
        assert!(
            (obs.range().value() - expected_range_km).abs() < 1e-3,
            "range={}",
            obs.range()
        );
        let _ = KM;
    }

    #[test]
    fn azimuth_singularity_is_scale_free_and_preserves_near_zenith() {
        // Exact vertical (any scale) → az = 0.
        assert_eq!(azimuth_from_enu(0.0, 0.0, 1.0), 0.0);
        assert_eq!(azimuth_from_enu(0.0, 0.0, -1e6), 0.0);
        // Noise at/below the ENU orthonormality floor is singular.
        assert_eq!(azimuth_from_enu(1e-13, -1e-13, 1.0), 0.0);
        assert_eq!(azimuth_from_enu(4.0 * f64::EPSILON, 0.0, 1.0), 0.0);
        // A distinguishable near-zenith offset keeps its quadrant.
        // horiz ≈ 1e-8 ≫ 1e-12; expected az = atan2(e,n) = 90°.
        let az = azimuth_from_enu(1e-8, 0.0, 1.0);
        assert!((az - core::f64::consts::FRAC_PI_2).abs() < 1e-12, "az={az}");
        let (alt, az) = angles_from_horizontal_cartesian(0.0, 1e-8, 1.0);
        assert!(alt.value() > 89.999, "alt={}", alt.value());
        assert!((az.value() - 90.0).abs() < 1e-6, "az={}", az.value());
    }

    #[test]
    fn direction_roundtrip_lunar() {
        let site = SurfaceSite::<Moon>::new(0.0 * DEG, 0.0 * DEG, 0.0 * M).unwrap();
        let jd = crate::J2000;
        let ctx =
            AstroContext::<crate::coordinates::transform::context::DefaultEphemeris>::default();
        // Local east in body-fixed at (0,0): east = (0,1,0)
        let east_fixed = Direction::<MoonIauMeanEarth>::new_unchecked(0.0, 1.0, 0.0);
        let icrs: Direction<ICRS> = east_fixed.to_frame_with(&jd, &ctx);
        let horiz = direction_to_horizontal(&icrs, &site, &jd, &ctx).unwrap();
        let sph = horiz.to_spherical();
        assert!((sph.alt().value()).abs() < 1e-6, "alt={}", sph.alt());
        let az = sph.az().value();
        assert!((az - 90.0).abs() < 1e-4, "az={az}");

        let back: Direction<ICRS> = horizontal_to_direction(&horiz, &site, &jd, &ctx).unwrap();
        // acos(dot) loses all angular resolution near zero (O(sqrt(eps))).
        // The chord residual is O(angle) and resolves machine-scale errors.
        let residual = norm([
            icrs.x() - back.x(),
            icrs.y() - back.y(),
            icrs.z() - back.z(),
        ]);
        assert!(residual < 1e-10, "roundtrip chord residual={residual}");
    }
}
