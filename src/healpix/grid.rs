// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::spherical;
use crate::healpix::{HealpixError, HealpixIndex, HealpixOrdering, Nside, Result};
use std::f64::consts::PI;

use super::ring;

/// HEALPix grid definition.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct HealpixGrid {
    nside: Nside,
    ordering: HealpixOrdering,
}

impl HealpixGrid {
    /// Construct a validated HEALPix grid.
    pub fn new(nside: Nside, ordering: HealpixOrdering) -> Result<Self> {
        if matches!(ordering, HealpixOrdering::Nested) && !nside.get().is_power_of_two() {
            return Err(HealpixError::NestedNsideNotPowerOfTwo(nside.get()));
        }

        Ok(Self { nside, ordering })
    }

    /// Construct a RING-ordered grid.
    pub fn ring(nside: Nside) -> Result<Self> {
        Self::new(nside, HealpixOrdering::Ring)
    }

    /// Return the grid `nside`.
    #[must_use]
    pub fn nside(&self) -> Nside {
        self.nside
    }

    /// Return the grid ordering.
    #[must_use]
    pub fn ordering(&self) -> HealpixOrdering {
        self.ordering
    }

    /// Return the number of pixels, `12 * nside^2`.
    #[must_use]
    pub fn npix(&self) -> u64 {
        let nside = u64::from(self.nside.get());
        12 * nside * nside
    }

    /// Return the solid angle of each pixel in steradians.
    #[must_use]
    pub fn pixel_area_sr(&self) -> f64 {
        4.0 * PI / self.npix() as f64
    }

    /// Return the solid angle of each pixel in square degrees.
    #[must_use]
    pub fn pixel_area_deg2(&self) -> f64 {
        self.pixel_area_sr() * (180.0 / PI).powi(2)
    }

    /// Validate that `index` belongs to this grid.
    pub fn validate_index(&self, index: HealpixIndex) -> Result<()> {
        let npix = self.npix();
        if index.get() >= npix {
            return Err(HealpixError::PixelIndexOutOfRange {
                index: index.get(),
                npix,
            });
        }
        Ok(())
    }

    /// Convert a typed cartesian direction to a HEALPix pixel index.
    pub fn direction_to_pixel<F>(&self, direction: Direction<F>) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        let [x, y, z] = direction.as_array();
        let r_xy = x.hypot(y);
        let r = r_xy.hypot(z);
        if !r.is_finite() || r == 0.0 {
            return Err(HealpixError::InvalidAngles {
                reason: "direction vector must be finite and non-zero",
            });
        }
        let lon = y.atan2(x);
        let lat = (z / r).clamp(-1.0, 1.0).asin();
        self.lon_lat_rad_to_pixel_checked(lon, lat)
    }

    /// Convert a typed spherical direction to a HEALPix pixel index.
    ///
    /// Siderust spherical directions store latitude-like polar angle first and
    /// longitude-like azimuth second, both as typed angular quantities.
    pub fn spherical_to_pixel<F>(
        &self,
        direction: spherical::Direction<F>,
    ) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        self.lon_lat_rad_to_pixel_checked(
            direction.azimuth.value().to_radians(),
            direction.polar.value().to_radians(),
        )
    }

    /// Return the typed cartesian direction at the centre of a pixel.
    pub fn pixel_center<F>(&self, index: HealpixIndex) -> Result<Direction<F>>
    where
        F: ReferenceFrame,
    {
        let (lon, lat) = self.pixel_to_lon_lat_rad_checked(index)?;
        Ok(ring::direction_from_lon_lat_rad(lon, lat))
    }

    /// Crate-private raw RING kernel from longitude/latitude in radians.
    pub(crate) fn lon_lat_rad_to_pixel_checked(
        &self,
        lon_rad: f64,
        lat_rad: f64,
    ) -> Result<HealpixIndex> {
        ring::lon_lat_rad_to_pixel_checked(self, lon_rad, lat_rad)
    }

    /// Crate-private raw RING kernel to longitude/latitude in radians.
    pub(crate) fn pixel_to_lon_lat_rad_checked(
        &self,
        index: HealpixIndex,
    ) -> Result<(f64, f64)> {
        ring::pixel_to_lon_lat_rad_checked(self, index)
    }
}
