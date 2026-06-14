// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Local east-north-up (ENU) reference frame at an observer site.

use affn::cartesian::Direction;
use affn::frames::ReferenceFrame;
use qtty::angular::Radians;

/// Local east-north-up (ENU) basis at an observer site whose geodetic
/// latitude/longitude is known.
///
/// `lat` and `lon` are the geodetic latitude and longitude of the site
/// expressed in the same frame the observer/target positions are given in
/// (typically ECEF).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalFrame {
    /// Geodetic latitude.
    pub lat: Radians,
    /// Geodetic longitude.
    pub lon: Radians,
}

impl LocalFrame {
    /// Build a local frame from typed geodetic coordinates.
    pub const fn new(lat: Radians, lon: Radians) -> Self {
        Self { lat, lon }
    }

    /// Compute the east, north, and up unit-direction vectors for this local frame.
    pub(in crate::mission::geometry) fn basis<F: ReferenceFrame>(
        &self,
    ) -> (Direction<F>, Direction<F>, Direction<F>) {
        let (sp, cp) = self.lat.sin_cos();
        let (sl, cl) = self.lon.sin_cos();
        let east = Direction::new_unchecked(-sl, cl, 0.0);
        let north = Direction::new_unchecked(-sp * cl, -sp * sl, cp);
        let up = Direction::new_unchecked(cp * cl, cp * sl, sp);
        (east, north, up)
    }
}
