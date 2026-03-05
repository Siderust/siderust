// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Center Transformations
//!
//! This module provides transformations between different astronomical reference centers
//! for **position** types only.
//!
//! ## Mathematical Foundations
//!
//! Center transformations are translations in affine space - they change the origin
//! from which positions are measured. Only affine objects (positions) can undergo
//! center transformations.
//!
//! **Directions and velocities are free vectors and cannot be center-transformed.**
//! Attempting to "change the center" of a direction is mathematically undefined.
//!
//! ## For Observer-Dependent Directions
//!
//! To compute the direction to a target as seen from an observer, use
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) with two positions:
//!
//! ```rust
//! use siderust::coordinates::cartesian::{line_of_sight, Position};
//! use siderust::coordinates::centers::Geocentric;
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use qtty::*;
//!
//! let observer =
//!     Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(0.0, 0.0, 0.0);
//! let target =
//!     Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(1.0, 1.0, 1.0);
//!
//! let direction = line_of_sight(&observer, &target);
//! ```

pub mod position;

// Re-export expert API for topocentric transforms with custom context.
pub use position::to_topocentric::to_topocentric_with_ctx;

use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::*;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::CenterShiftProvider;
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, LengthUnit, Quantity};

// =============================================================================
// IntoTransformArgs — converts (params, jd) or just jd into the full argument
// =============================================================================

/// Converts a caller-supplied value into the `(params, jd)` pair required by
/// [`TransformCenter::to_center`].
///
/// This allows calling `to_center` with different argument styles:
///
/// | Call site | `args` type | Condition |
/// |---|---|---|
/// | `pos.to_center(jd)` | [`JulianDate`] | `C2::Params = ()` |
/// | `pos.to_center((site, jd))` | `(Geodetic<ECEF>, JulianDate)` | `C2 = Topocentric` |
/// | `pos.to_center((orbit_params, jd))` | `(BodycentricParams, JulianDate)` | `C2 = Bodycentric` |
pub trait IntoTransformArgs<Params> {
    /// Decompose into parameters and Julian date.
    fn into_params_jd(self) -> (Params, JulianDate);
}

/// For standard centers (`Params = ()`): passing just a [`JulianDate`] is enough.
impl IntoTransformArgs<()> for JulianDate {
    #[inline]
    fn into_params_jd(self) -> ((), JulianDate) {
        ((), self)
    }
}

/// For **any** center: passing a `(Params, JulianDate)` tuple always works.
///
/// This is the canonical form for non-trivial parameter types like
/// [`BodycentricParams`](crate::coordinates::centers::BodycentricParams) or
/// [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic).
impl<P: Clone> IntoTransformArgs<P> for (P, JulianDate) {
    #[inline]
    fn into_params_jd(self) -> (P, JulianDate) {
        self
    }
}

/// Trait for transforming a [`Position`] to a different reference center.
///
/// This is the **single, canonical API** for all center shifts.
///
/// | Target center | Call site | Notes |
/// |---|---|---|
/// | Barycentric / Heliocentric / Geocentric (identity too) | `pos.to_center(jd)` | Pass only the [`JulianDate`] — no `()` |
/// | Bodycentric | `pos.to_center((orbit_params, jd))` | [`BodycentricParams`](crate::coordinates::centers::BodycentricParams) |
/// | Topocentric | `pos.to_center((site, jd))` | [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic) |
/// | Bodycentric → Geocentric | `bary.to_center(jd)` | Same as standard |
/// | Topocentric → Geocentric | `topo.to_center(jd)` | Same as standard |
///
/// # Type Parameters
///
/// - `C2`: Target reference center.  Its [`ReferenceCenter::Params`] type
///   determines what `args` must be passed.
/// - `F`: Shared reference frame (center shifts are frame-respecting translations).
/// - `U`: Length unit.
pub trait TransformCenter<C2: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Transform to the target center using the default ephemeris / EOP.
    ///
    /// `args` is either a bare [`JulianDate`] (for standard centers whose
    /// `Params = ()`) or a `(params, jd)` tuple for parameterised centers.
    /// See [`IntoTransformArgs`] for all supported forms.
    fn to_center<A: IntoTransformArgs<C2::Params>>(&self, args: A) -> Position<C2, F, U> {
        let (params, jd) = args.into_params_jd();
        self.to_center_with(params, jd, &AstroContext::default())
    }

    /// Transform to the target center with a custom [`AstroContext`].
    ///
    /// Use this to override the ephemeris, EOP, or nutation model.
    fn to_center_with(
        &self,
        params: C2::Params,
        jd: JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F, U>;
}

// =============================================================================
// Blanket impl for standard (Params = ()) centers via CenterShiftProvider
// =============================================================================

/// Blanket implementation covering all standard-center pairs
/// (Barycentric ↔ Heliocentric ↔ Geocentric, and identity).
///
/// The shift vector is looked up from [`CenterShiftProvider`], which
/// consults the configured ephemeris. No extra concrete per-pair impl is
/// needed; deleting the old `to_barycentric.rs`, `to_geocentric.rs`, and
/// `to_heliocentric.rs` files is intentional.
impl<C1, C2, F, U> TransformCenter<C2, F, U> for Position<C1, F, U>
where
    C1: ReferenceCenter<Params = ()>,
    C2: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    U: LengthUnit,
    (): CenterShiftProvider<C1, C2, F>,
{
    fn to_center_with(
        &self,
        _params: (),
        jd: JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F, U> {
        let shift = <() as CenterShiftProvider<C1, C2, F>>::shift(jd, ctx);
        let sx = Quantity::<AstronomicalUnit>::new(shift[0]).to::<U>();
        let sy = Quantity::<AstronomicalUnit>::new(shift[1]).to::<U>();
        let sz = Quantity::<AstronomicalUnit>::new(shift[2]).to::<U>();
        Position::new(self.x() + sx, self.y() + sy, self.z() + sz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calculus::ephemeris::Ephemeris;
    use crate::coordinates::cartesian;
    use crate::coordinates::transform::context::DefaultEphemeris;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use crate::time::JulianDate;
    use qtty::AstronomicalUnit;

    const EPSILON: f64 = 1e-8;

    #[test]
    fn test_position_barycentric_to_geocentric() {
        let earth_bary = DefaultEphemeris::earth_barycentric(JulianDate::J2000);
        let earth_geo: cartesian::position::EclipticMeanJ2000<AstronomicalUnit, Geocentric> =
            earth_bary.transform(JulianDate::J2000);
        let expected_earth_geo =
            cartesian::position::EclipticMeanJ2000::<AstronomicalUnit, Geocentric>::CENTER;
        assert_cartesian_eq!(
            &earth_geo,
            &expected_earth_geo,
            EPSILON,
            "Earth in Geocentric should be at origin. Current: {:?}",
            earth_geo
        );
    }

    #[test]
    fn test_position_heliocentric_to_geocentric() {
        let earth_helio = DefaultEphemeris::earth_heliocentric(JulianDate::J2000);
        let earth_geo: cartesian::position::EclipticMeanJ2000<AstronomicalUnit, Geocentric> =
            earth_helio.transform(JulianDate::J2000);
        let expected =
            cartesian::position::EclipticMeanJ2000::<AstronomicalUnit, Geocentric>::CENTER;
        assert_cartesian_eq!(&earth_geo, &expected, EPSILON);
    }
}
