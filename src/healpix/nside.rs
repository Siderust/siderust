// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! HEALPix `nside` resolution parameter.
//!
//! HEALPix divides the sphere into `12 * nside^2` equal-area pixels. `nside`
//! must be strictly positive; NESTED indexing additionally requires a
//! power-of-two value and is validated by [`crate::healpix::HealpixGrid`].

use crate::healpix::{HealpixError, Result};

/// Largest `nside` whose `12 * nside^2` pixel count fits in `u64`.
pub const MAX_NSIDE: u32 = 1_239_850_262;

/// Validated HEALPix `nside` resolution parameter.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Nside(u32);

impl Nside {
    /// Construct a validated `nside`.
    pub fn new(value: u32) -> Result<Self> {
        if value == 0 {
            return Err(HealpixError::InvalidNside);
        }
        if value > MAX_NSIDE {
            return Err(HealpixError::NsideTooLarge(value));
        }
        Ok(Self(value))
    }

    /// Return the underlying integer value.
    #[must_use]
    pub fn get(self) -> u32 {
        self.0
    }
}

impl TryFrom<u32> for Nside {
    type Error = HealpixError;

    fn try_from(value: u32) -> Result<Self> {
        Self::new(value)
    }
}

impl From<Nside> for u32 {
    fn from(value: Nside) -> Self {
        value.0
    }
}
