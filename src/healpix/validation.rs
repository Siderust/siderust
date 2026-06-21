// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixMap, Result};

/// Validate that a HEALPix map has one value for every pixel.
pub fn validate_healpix_map_complete<F, T>(map: &HealpixMap<F, T>) -> Result<()>
where
    F: ReferenceFrame,
{
    let len = map.values().len();
    let npix = map.grid().npix();
    if u64::try_from(len).expect("usize length fits u64") != npix {
        return Err(HealpixError::MapLengthMismatch { len, npix });
    }
    Ok(())
}
