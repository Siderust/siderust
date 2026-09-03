// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! HEALPix grid definition and frame-neutral pixelization operations.
//!
//! A [`HealpixGrid`] stores the HEALPix resolution parameter `nside` and the
//! selected ordering scheme. The grid itself is intentionally not tied to an
//! astronomical frame; frame semantics are carried by typed input and output
//! directions.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::spherical;
use crate::healpix::{HealpixError, HealpixIndex, HealpixOrdering, Nside, Result};
use std::f64::consts::PI;

use super::ring;

/// HEALPix grid descriptor.
///
/// The grid is a small value type containing only the validated `nside`
/// resolution and the ordering scheme. It does not store pixel values and can
/// therefore be copied cheaply. Pixel values are owned by [`HealpixMap`](crate::healpix::HealpixMap).
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct HealpixGrid {
    nside: Nside,
    ordering: HealpixOrdering,
}

impl HealpixGrid {
    /// Construct a validated HEALPix grid.
    ///
    /// NESTED ordering currently requires a power-of-two `nside`, matching the
    /// hierarchical subdivision used by the NESTED indexing scheme.
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

    /// Return the grid resolution parameter.
    #[must_use]
    pub fn nside(&self) -> Nside {
        self.nside
    }

    /// Return the pixel ordering scheme.
    #[must_use]
    pub fn ordering(&self) -> HealpixOrdering {
        self.ordering
    }

    /// Return the number of pixels in the grid, `12 * nside^2`.
    #[must_use]
    pub fn npix(&self) -> u64 {
        let nside = u64::from(self.nside.get());
        12 * nside * nside
    }

    /// Return the equal-area solid angle of each pixel in steradians.
    #[must_use]
    pub fn pixel_area_sr(&self) -> f64 {
        4.0 * PI / self.npix() as f64
    }

    /// Return the equal-area solid angle of each pixel in square degrees.
    #[must_use]
    pub fn pixel_area_deg2(&self) -> f64 {
        self.pixel_area_sr() * (180.0 / PI).powi(2)
    }

    /// Validate that a pixel index belongs to this grid.
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

    /// Convert a typed direction to a HEALPix pixel index.
    ///
    /// The direction frame is preserved at the type level. The implementation
    /// treats the direction as a unit vector in the axes of that frame and never
    /// exposes frame-specific angular names such as right ascension,
    /// declination, longitude, or latitude.
    pub fn direction_to_pixel<F>(&self, direction: Direction<F>) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        ring::unit_vector_to_pixel(self, direction.as_array())
    }

    /// Return the typed direction at the center of a pixel.
    ///
    /// The returned direction is expressed in the requested reference frame type
    /// `F`; callers must choose a frame consistent with the map or grid they are
    /// working with.
    pub fn pixel_center<F>(&self, index: HealpixIndex) -> Result<Direction<F>>
    where
        F: ReferenceFrame,
    {
        let (theta, phi) = ring::pixel_to_theta_phi(self, index)?;
        Ok(ring::direction_from_theta_phi(theta, phi))
    }

    /// Return the center of a pixel as a typed spherical direction.
    ///
    /// The returned [`spherical::Direction`] is expressed in the requested
    /// reference frame type `F`; the grid itself is frame-neutral, so callers
    /// must choose a frame consistent with the map or grid they are working
    /// with.
    ///
    /// HEALPix describes a center using `theta`, its colatitude, and `phi`, its
    /// longitude. This method maps those angles directly to Siderust's
    /// spherical convention: `polar = π/2 - theta` and `azimuth = phi`.
    /// The result uses the canonical angular ranges supplied by `affn` and
    /// `qtty`: polar is in `[-90°, +90°]` and azimuth is in `[0°, 360°)`.
    ///
    /// No Cartesian conversion is performed.
    pub fn pixel_center_spherical<F>(&self, index: HealpixIndex) -> Result<spherical::Direction<F>>
    where
        F: ReferenceFrame,
    {
        let (theta, phi) = ring::pixel_to_theta_phi(self, index)?;
        Ok(ring::spherical_direction_from_theta_phi(theta, phi))
    }
}
