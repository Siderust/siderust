// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Local east-north-up (ENU) reference frame at an observer site.

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

    /// Compute the east, north, and up unit vectors for this local frame.
    pub(in crate::mission::geometry) fn basis(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (sp, cp) = (self.lat.value().sin(), self.lat.value().cos());
        let (sl, cl) = (self.lon.value().sin(), self.lon.value().cos());
        let east = [-sl, cl, 0.0];
        let north = [-sp * cl, -sp * sl, cp];
        let up = [cp * cl, cp * sl, sp];
        (east, north, up)
    }
}
