// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::healpix::HealpixOrdering;

/// Result alias for HEALPix operations.
pub type Result<T> = std::result::Result<T, HealpixError>;

/// Error type for HEALPix grid, index, and map operations.
#[derive(Debug, thiserror::Error, Clone, PartialEq)]
pub enum HealpixError {
    /// `nside` must be greater than zero.
    #[error("HEALPix nside must be greater than zero")]
    InvalidNside,
    /// `nside` is too large for safe pixel-count arithmetic.
    #[error("HEALPix nside {0} is too large")]
    NsideTooLarge(u32),
    /// NESTED ordering requires a power-of-two `nside`.
    #[error("HEALPix NESTED ordering requires power-of-two nside, got {0}")]
    NestedNsideNotPowerOfTwo(u32),
    /// The selected ordering is not yet implemented for this operation.
    #[error("HEALPix ordering {0:?} is not supported by this operation yet")]
    UnsupportedOrdering(HealpixOrdering),
    /// A pixel index was outside the grid range.
    #[error("HEALPix pixel index {index} is outside the valid range 0..{npix}")]
    PixelIndexOutOfRange {
        /// Invalid pixel index.
        index: u64,
        /// Number of pixels in the grid.
        npix: u64,
    },
    /// The provided longitude/latitude was invalid.
    #[error("invalid HEALPix angular input: {reason}")]
    InvalidAngles {
        /// Explanation of the angular validation failure.
        reason: &'static str,
    },
    /// A map has a number of values inconsistent with its grid.
    #[error("HEALPix map length {len} does not match grid npix {npix}")]
    MapLengthMismatch {
        /// Number of stored values.
        len: usize,
        /// Number of pixels expected by the grid.
        npix: u64,
    },
    /// A finite, non-negative map value check failed.
    #[error("HEALPix map contains an invalid non-finite or negative value")]
    InvalidMapValue,
}
