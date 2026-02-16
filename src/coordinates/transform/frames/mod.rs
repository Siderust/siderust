// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

mod bias;
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
use qtty::LengthUnit;

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
        Position::from_vec3(self.center_params().clone(), *self.as_vec3())
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
