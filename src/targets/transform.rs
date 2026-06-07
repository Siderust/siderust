// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Coordinate-type conversions for [`CoordinateWithPM`]
//!
//! ## Scientific scope
//!
//! Astrometric pipelines frequently need to express the same catalog position
//! in different reference frames or reference centers: e.g. transforming from
//! a heliocentric J2000 ecliptic catalog position to a geocentric equatorial
//! of-date position for comparison with an observation. This module provides
//! blanket [`From`] implementations so those transformations compose
//! transparently with the `CoordinateWithPM<T>` container, preserving the
//! epoch and proper-motion model across the conversion.
//!
//! ## Technical scope
//!
//! Two blanket impls are provided:
//!
//! - `From<&CoordinateWithPM<cartesian::Position<C1,F1,U>>>`
//!   `for CoordinateWithPM<cartesian::Position<C2,F2,U>>` — chains a frame
//!   transform and a center transform, both evaluated at the stored epoch.
//!
//! - `From<&CoordinateWithPM<spherical::Position<C1,F1,U>>>`
//!   `for CoordinateWithPM<spherical::Position<C2,F2,U>>` — converts to
//!   Cartesian, applies both transforms, then converts back to spherical.
//!
//! Both impls delegate to the [`Transform`](crate::coordinates::transform::Transform)
//! trait and are gated on the same bounds, so the compiler statically checks
//! that the requested center/frame combination has a valid transform path.
//!
//! ## References
//!
//! - IAU SOFA Library, function `iauPmpx` (proper-motion and parallax).
//! - Seidelmann, P. K. (1992). *Explanatory Supplement to the Astronomical
//!   Almanac*, §3.3. University Science Books.

use super::CoordinateWithPM;
use crate::coordinates::{cartesian, centers::*, frames::*, spherical, transform::Transform};
use crate::qtty::LengthUnit;

/// Blanket implementation to allow chaining two consecutive `Transform` operations.
///
/// This implementation allows converting a `CoordinateWithPM` in Cartesian coordinates from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U> From<&CoordinateWithPM<cartesian::Position<C1, F1, U>>>
    for CoordinateWithPM<cartesian::Position<C2, F2, U>>
where
    cartesian::Position<C1, F1, U>: Transform<cartesian::Position<C1, F2, U>>, // transform frame
    cartesian::Position<C1, F2, U>: Transform<cartesian::Position<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &CoordinateWithPM<cartesian::Position<C1, F1, U>>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        // Step 2: Transform to new center, now using the new frame.
        Self::new_raw(
            orig.position.transform(orig.time).transform(orig.time),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}

/// Blanket implementation for transforming `CoordinateWithPM` in spherical coordinates,
/// involving frame and center changes. Internally uses Cartesian conversions.
///
/// The transformation follows these steps:
/// 1. Convert spherical coordinates to Cartesian.
/// 2. Apply frame transformation.
/// 3. Apply center transformation.
/// 4. Convert back to spherical coordinates.
impl<C1, F1, C2, F2, U> From<&CoordinateWithPM<spherical::Position<C1, F1, U>>>
    for CoordinateWithPM<spherical::Position<C2, F2, U>>
where
    cartesian::Position<C1, F1, U>: Transform<cartesian::Position<C1, F2, U>>, // transform frame
    cartesian::Position<C1, F2, U>: Transform<cartesian::Position<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    U: LengthUnit,
{
    fn from(orig: &CoordinateWithPM<spherical::Position<C1, F1, U>>) -> Self {
        // Step 1: Convert spherical to Cartesian
        // Step 2: Transform to new frame
        // Step 3: Transform to new center
        // Step 4: Convert back to spherical
        Self::new_raw(
            spherical::Position::from_cartesian(
                &orig
                    .position
                    .to_cartesian()
                    .transform(orig.time)
                    .transform(orig.time),
            ),
            orig.time,
            orig.proper_motion.clone(),
        )
    }
}
