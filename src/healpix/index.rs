// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! HEALPix pixel index type.
//!
//! Pixel indices are zero-based offsets into the ordering scheme selected by a
//! [`crate::healpix::HealpixGrid`]. The index type stores the integer value,
//! while grid-range validation is performed against a concrete grid.

/// HEALPix pixel index.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct HealpixIndex(u64);

impl HealpixIndex {
    /// Construct a pixel index without grid-range validation.
    ///
    /// Use [`crate::healpix::HealpixGrid::validate_index`] when a concrete grid is available.
    #[must_use]
    pub fn new(value: u64) -> Self {
        Self(value)
    }

    /// Return the underlying integer index.
    #[must_use]
    pub fn get(self) -> u64 {
        self.0
    }
}

impl From<HealpixIndex> for u64 {
    fn from(value: HealpixIndex) -> Self {
        value.0
    }
}
