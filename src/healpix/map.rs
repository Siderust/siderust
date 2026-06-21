//! Frame-typed HEALPix map storage.
//!
//! A HEALPix map is represented as one value per pixel in the ordering defined
//! by its [`HealpixGrid`]. The coordinate frame is carried by the type parameter
//! `F`, while the stored value type `T` remains domain-specific.

use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixGrid, Result};
use std::marker::PhantomData;

/// Complete HEALPix map over a specific reference frame.
///
/// The constructor enforces that the number of stored values matches
/// `grid.npix()`. After construction, callers can mutate values in place but
/// cannot resize the map through the public API.
#[derive(Debug, Clone, PartialEq)]
pub struct HealpixMap<F, T>
where
    F: ReferenceFrame,
{
    grid: HealpixGrid,
    values: Vec<T>,
    marker: PhantomData<F>,
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
            marker: PhantomData,
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
    ///
    /// The returned slice cannot change the map length, so the map-completeness
    /// invariant is preserved.
    #[must_use]
    pub fn values_mut(&mut self) -> &mut [T] {
        &mut self.values
    }

    /// Consume the map and return its values.
    #[must_use]
    pub fn into_values(self) -> Vec<T> {
        self.values
    }
}
