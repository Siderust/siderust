// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixGrid, Result};
use std::marker::PhantomData;

/// Frame-typed HEALPix map container.
#[derive(Debug, Clone, PartialEq)]
pub struct HealpixMap<F, T>
where
    F: ReferenceFrame,
{
    grid: HealpixGrid,
    values: Vec<T>,
    frame: PhantomData<F>,
}

impl<F, T> HealpixMap<F, T>
where
    F: ReferenceFrame,
{
    /// Construct a map and validate that `values.len() == grid.npix()`.
    pub fn new(grid: HealpixGrid, values: Vec<T>) -> Result<Self> {
        let len = values.len();
        let npix = grid.npix();
        if u64::try_from(len).expect("usize length fits u64") != npix {
            return Err(HealpixError::MapLengthMismatch { len, npix });
        }
        Ok(Self {
            grid,
            values,
            frame: PhantomData,
        })
    }

    /// Return the grid shared by all pixel values.
    #[must_use]
    pub fn grid(&self) -> HealpixGrid {
        self.grid
    }

    /// Return all pixel values in grid ordering.
    #[must_use]
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// Return all pixel values in grid ordering for in-place value updates.
    #[must_use]
    pub fn values_mut(&mut self) -> &mut [T] {
        &mut self.values
    }

    /// Consume the map and return its values.
    #[must_use]
    pub fn into_values(self) -> Vec<T> {
        self.values
    }

    /// Return the frame marker carried by this map type.
    #[must_use]
    pub fn frame_marker(&self) -> PhantomData<F> {
        self.frame
    }
}
