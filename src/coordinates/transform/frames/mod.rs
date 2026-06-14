// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame rotation trait and blanket implementations.
//!
//! ## Scientific scope
//!
//! Frame rotations change the orientation of coordinate axes without moving the
//! origin. The reference frames bridged in this sub-tree include:
//!
//! - **ICRS / GCRS** — pseudo-inertial, tied to extragalactic radio sources.
//! - **EquatorialMeanJ2000 / EclipticMeanJ2000** — J2000.0 epoch realisation.
//! - **EquatorialMeanOfDate / EquatorialTrueOfDate** — precession and nutation.
//! - **CIRS / TIRS / ITRF / ECEF** — Earth-fixed chain (ERA, polar motion).
//! - **Horizontal** — local horizon frame, tied to an observer site.
//!
//! ## Technical scope
//!
//! [`TransformFrame`] is the low-level, statically-dispatched interface for
//! frame rotations. Each `impl` in the child modules applies a fixed or
//! time-dependent 3×3 rotation matrix to the underlying Cartesian vector.
//!
//! Blanket impls here route spherical types through a Cartesian round-trip so
//! that each frame-pair only needs a single Cartesian `impl`.
//!
//! The higher-level, runtime API — which selects the nutation model at
//! run-time and carries an `AstroContext` — lives in
//! [`crate::coordinates::transform::ext`] and uses `FrameRotationProvider`.
//!
//! ## References
//!
//! - IAU SOFA Tools for Earth Attitude (2023): <https://www.iausofa.org>
//! - IERS Conventions 2010 (IERS Technical Note No. 36).

#![allow(unreachable_pub, missing_docs)]

pub(crate) mod bias;
pub mod direction;
pub mod to_ecliptic;
pub mod to_equatorial;
pub mod to_horizontal;
pub mod to_icrs;
pub mod velocity;

use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::{cartesian, spherical};
use crate::qtty::LengthUnit;

pub trait TransformFrame<Coord> {
    fn to_frame(&self) -> Coord;
}

// Implement Identity frame transform
impl<C, F, U> TransformFrame<Position<C, F, U>> for Position<C, F, U>
where
    U: LengthUnit,
    F: MutableFrame,
    C: ReferenceCenter,
{
    fn to_frame(&self) -> Position<C, F, U> {
        Position::from_array(self.center_params().clone(), *self.as_array())
    }
}

impl<C, F1, F2, U> TransformFrame<spherical::Position<C, F2, U>> for spherical::Position<C, F1, U>
where
    cartesian::Position<C, F1, U>: TransformFrame<cartesian::Position<C, F2, U>>,
    C: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn to_frame(&self) -> spherical::Position<C, F2, U> {
        spherical::Position::from_cartesian(&self.to_cartesian().to_frame())
    }
}

/// Blanket `TransformFrame` for spherical **directions** via cartesian round-trip.
///
/// Any pair `(F1, F2)` for which `cartesian::Direction<F1>: TransformFrame<cartesian::Direction<F2>>`
/// is satisfied automatically gets a spherical direction transform.
impl<F1, F2> TransformFrame<spherical::Direction<F2>> for spherical::Direction<F1>
where
    cartesian::Direction<F1>: TransformFrame<cartesian::Direction<F2>>,
    F1: MutableFrame,
    F2: MutableFrame,
{
    fn to_frame(&self) -> spherical::Direction<F2> {
        spherical::Direction::from_cartesian(&self.to_cartesian().to_frame())
    }
}

// Note: The to_frame() method for spherical::Position and spherical::Direction
// is provided by the TransformFrame trait implementations above. Users can call:
//   sph_pos.to_frame::<NewFrame>()
//   sph_dir.to_frame::<NewFrame>()
