// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generated stellar surface-brightness map wrapper.
//!
//! The wrapper keeps the typed Galactic HEALPix map and its provenance metadata
//! together as one generated data product.

use crate::coordinates::frames::Galactic;
use crate::healpix::{HealpixGrid, HealpixMap};
use crate::starlight::{StellarMapProvenance, StellarSurfaceBrightness};

/// Generated Galactic stellar surface-brightness map with provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct StellarSurfaceBrightnessMap {
    map: HealpixMap<Galactic, StellarSurfaceBrightness>,
    provenance: StellarMapProvenance,
}

impl StellarSurfaceBrightnessMap {
    /// Construct a generated stellar map from validated map data and provenance.
    #[must_use]
    pub fn new(
        map: HealpixMap<Galactic, StellarSurfaceBrightness>,
        provenance: StellarMapProvenance,
    ) -> Self {
        Self { map, provenance }
    }

    /// Return the underlying HEALPix grid.
    #[must_use]
    pub fn grid(&self) -> HealpixGrid {
        self.map.grid()
    }

    /// Return the pixel values in grid ordering.
    #[must_use]
    pub fn values(&self) -> &[StellarSurfaceBrightness] {
        self.map.values()
    }

    /// Return the provenance metadata attached to this generated map.
    #[must_use]
    pub fn provenance(&self) -> &StellarMapProvenance {
        &self.provenance
    }

    /// Return the inner typed HEALPix map.
    #[must_use]
    pub fn healpix_map(&self) -> &HealpixMap<Galactic, StellarSurfaceBrightness> {
        &self.map
    }

    /// Consume the wrapper and return the map and provenance.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        HealpixMap<Galactic, StellarSurfaceBrightness>,
        StellarMapProvenance,
    ) {
        (self.map, self.provenance)
    }
}
