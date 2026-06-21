// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::Galactic;
use crate::coordinates::transform::TransformFrame;
use crate::healpix::{HealpixGrid, HealpixMap};
use crate::starlight::{
    flux_10mag_units, validate_stellar_map_values, ApparentMagnitude, Result,
    StellarCatalogueRecord, StellarMapError, StellarMapProvenance, StellarSurfaceBrightness,
};

/// Builder for Galactic stellar surface-brightness HEALPix maps.
#[derive(Debug, Clone, PartialEq)]
pub struct StellarSurfaceBrightnessMapBuilder {
    /// Output HEALPix grid.
    pub grid: HealpixGrid,
    /// Optional inclusive lower V-magnitude cut.
    pub min_v_mag: Option<ApparentMagnitude>,
    /// Optional inclusive upper V-magnitude cut.
    pub max_v_mag: Option<ApparentMagnitude>,
    /// Scale from V S10 per square degree to integrated photon radiance.
    pub integrated_per_v_s10: f64,
    /// Optional smoothing FWHM in degrees. The v1 builder requires `None`.
    pub smoothing_fwhm_deg: Option<f64>,
}

impl StellarSurfaceBrightnessMapBuilder {
    /// Build a Galactic HEALPix stellar surface-brightness map.
    pub fn build<I>(
        &self,
        records: I,
        _provenance: StellarMapProvenance,
    ) -> Result<HealpixMap<Galactic, StellarSurfaceBrightness>>
    where
        I: IntoIterator<Item = StellarCatalogueRecord>,
    {
        if !self.integrated_per_v_s10.is_finite() || self.integrated_per_v_s10 < 0.0 {
            return Err(StellarMapError::InvalidRadianceScale(
                self.integrated_per_v_s10,
            ));
        }
        if self.smoothing_fwhm_deg.is_some() {
            return Err(StellarMapError::UnsupportedSmoothing);
        }

        let npix = usize::try_from(self.grid.npix()).expect("HEALPix npix fits usize");
        let mut values = vec![StellarSurfaceBrightness::zero(); npix];
        let mut included = 0_u64;
        let mut contributed_flux = false;

        for record in records {
            if !record.weight.is_finite() || record.weight < 0.0 {
                return Err(StellarMapError::InvalidWeight(record.weight));
            }
            if !self.passes_v_magnitude_cut(record.v_mag) {
                continue;
            }

            let galactic: Direction<Galactic> = record.direction.to_frame();
            let index = self.grid.direction_to_pixel(galactic)?;
            let slot = usize::try_from(index.get()).expect("pixel index fits usize");
            included += 1;

            if let Some(mag) = record.b_mag {
                let flux = flux_10mag_units(mag) * record.weight;
                values[slot].b_s10 += flux;
                contributed_flux |= flux > 0.0;
            }
            if let Some(mag) = record.v_mag {
                let flux = flux_10mag_units(mag) * record.weight;
                values[slot].v_s10 += flux;
                contributed_flux |= flux > 0.0;
            }
        }

        if included == 0 {
            return Err(StellarMapError::EmptyFilteredCatalogue);
        }
        if !contributed_flux {
            return Err(StellarMapError::NoUsablePhotometry);
        }

        let pixel_area_deg2 = self.grid.pixel_area_deg2();
        for value in &mut values {
            value.b_s10 /= pixel_area_deg2;
            value.v_s10 /= pixel_area_deg2;
            value.integrated_ph_cm2_ns_sr = value.v_s10 * self.integrated_per_v_s10;
        }

        let map = HealpixMap::<Galactic, StellarSurfaceBrightness>::new(self.grid, values)?;
        validate_stellar_map_values(&map)?;
        Ok(map)
    }

    fn passes_v_magnitude_cut(&self, magnitude: Option<ApparentMagnitude>) -> bool {
        match magnitude {
            Some(value) => {
                self.min_v_mag.is_none_or(|min| value.value() >= min.value())
                    && self.max_v_mag.is_none_or(|max| value.value() <= max.value())
            }
            None => self.min_v_mag.is_none() && self.max_v_mag.is_none(),
        }
    }
}
