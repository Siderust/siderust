use crate::coordinates::frames::Galactic;
use crate::healpix::{HealpixGrid, HealpixMap};
use crate::starlight::{StellarMapProvenance, StellarSurfaceBrightness};

#[derive(Debug, Clone, PartialEq)]
pub struct StellarSurfaceBrightnessMap {
    map: HealpixMap<Galactic, StellarSurfaceBrightness>,
    provenance: StellarMapProvenance,
}

impl StellarSurfaceBrightnessMap {
    #[must_use]
    pub fn new(
        map: HealpixMap<Galactic, StellarSurfaceBrightness>,
        provenance: StellarMapProvenance,
    ) -> Self {
        Self { map, provenance }
    }

    #[must_use]
    pub fn grid(&self) -> HealpixGrid {
        self.map.grid()
    }

    #[must_use]
    pub fn values(&self) -> &[StellarSurfaceBrightness] {
        self.map.values()
    }

    #[must_use]
    pub fn provenance(&self) -> &StellarMapProvenance {
        &self.provenance
    }

    #[must_use]
    pub fn healpix_map(&self) -> &HealpixMap<Galactic, StellarSurfaceBrightness> {
        &self.map
    }

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
