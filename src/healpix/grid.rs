use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixIndex, HealpixOrdering, Nside, Result};
use std::f64::consts::PI;

use super::ring;

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct HealpixGrid {
    nside: Nside,
    ordering: HealpixOrdering,
}

impl HealpixGrid {
    pub fn new(nside: Nside, ordering: HealpixOrdering) -> Result<Self> {
        if matches!(ordering, HealpixOrdering::Nested) && !nside.get().is_power_of_two() {
            return Err(HealpixError::NestedNsideNotPowerOfTwo(nside.get()));
        }
        Ok(Self { nside, ordering })
    }

    pub fn ring(nside: Nside) -> Result<Self> {
        Self::new(nside, HealpixOrdering::Ring)
    }

    #[must_use]
    pub fn nside(&self) -> Nside {
        self.nside
    }

    #[must_use]
    pub fn ordering(&self) -> HealpixOrdering {
        self.ordering
    }

    #[must_use]
    pub fn npix(&self) -> u64 {
        let nside = u64::from(self.nside.get());
        12 * nside * nside
    }

    #[must_use]
    pub fn pixel_area_sr(&self) -> f64 {
        4.0 * PI / self.npix() as f64
    }

    #[must_use]
    pub fn pixel_area_deg2(&self) -> f64 {
        self.pixel_area_sr() * (180.0 / PI).powi(2)
    }

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

    pub fn direction_to_pixel<F>(&self, direction: Direction<F>) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        ring::unit_vector_to_pixel(self, direction.as_array())
    }

    pub fn pixel_center<F>(&self, index: HealpixIndex) -> Result<Direction<F>>
    where
        F: ReferenceFrame,
    {
        let (theta, phi) = ring::pixel_to_theta_phi(self, index)?;
        Ok(ring::direction_from_theta_phi(theta, phi))
    }
}
