// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Stellar surface-brightness map construction.
//!
//! ## Scientific scope
//!
//! Integrated starlight is represented as a Galactic HEALPix map whose pixels
//! accumulate catalogue-star flux in tenth-magnitude-star units. The current
//! model is an explicitly labelled v1 approximation intended for reproducible
//! data-asset generation and validation, not full passband/spectral synthesis.
//!
//! ## Technical scope
//!
//! This module consumes in-memory catalogue records, transforms their typed
//! EquatorialMeanJ2000 directions to Galactic directions with Siderust frame
//! transforms, bins them into a typed [`HealpixMap<Galactic, StellarSurfaceBrightness>`],
//! and provides reusable validators for flux conservation and simple full-sky
//! sanity checks.
//!
//! ## References
//!
//! - Gorski, K. M. et al. (2005). *HEALPix: A Framework for High-Resolution
//!   Discretization and Fast Analysis of Data Distributed on the Sphere*, ApJ,
//!   622, 759.
//! - Leinert, C. et al. (1998). *The 1997 reference of diffuse night sky
//!   brightness*, A&AS, 127, 1.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{EquatorialMeanJ2000, Galactic};
use crate::coordinates::transform::TransformFrame;
use crate::healpix::{HealpixError, HealpixGrid, HealpixIndex, HealpixMap};

/// Result alias for stellar surface-brightness map operations.
pub type Result<T> = std::result::Result<T, StellarMapError>;

/// Error type for stellar surface-brightness map construction and validation.
#[derive(Debug, thiserror::Error)]
pub enum StellarMapError {
    /// Wrapped HEALPix error.
    #[error(transparent)]
    Healpix(#[from] HealpixError),
    /// A magnitude value was invalid.
    #[error("apparent magnitude must be finite, got {0}")]
    InvalidMagnitude(f64),
    /// A catalogue weight was invalid.
    #[error("catalogue weight must be finite and non-negative, got {0}")]
    InvalidWeight(f64),
    /// A radiance scale was invalid.
    #[error("integrated_per_v_s10 must be finite and non-negative, got {0}")]
    InvalidRadianceScale(f64),
    /// Smoothing was requested but is not implemented in the v1 builder.
    #[error("stellar-map smoothing is not implemented in the v1 builder")]
    UnsupportedSmoothing,
    /// No catalogue records survived filtering.
    #[error("no stellar catalogue records survived filtering")]
    EmptyFilteredCatalogue,
    /// A validation check failed.
    #[error("stellar map validation failed: {0}")]
    Validation(&'static str),
}

/// Apparent magnitude newtype used by the v1 photometric model.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct ApparentMagnitude(f64);

impl ApparentMagnitude {
    /// Construct a finite apparent magnitude.
    pub fn new(value: f64) -> Result<Self> {
        if !value.is_finite() {
            return Err(StellarMapError::InvalidMagnitude(value));
        }
        Ok(Self(value))
    }

    /// Return the magnitude value.
    #[must_use]
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Generic stellar catalogue record consumed by the map builder.
#[derive(Debug, Clone, PartialEq)]
pub struct StellarCatalogueRecord {
    /// Optional source identifier from the upstream catalogue.
    pub source_id: Option<String>,
    /// Source direction in the EquatorialMeanJ2000 frame.
    pub direction: Direction<EquatorialMeanJ2000>,
    /// Optional B-band apparent magnitude.
    pub b_mag: Option<ApparentMagnitude>,
    /// Optional V-band apparent magnitude.
    pub v_mag: Option<ApparentMagnitude>,
    /// Multiplicative catalogue weight; must be finite and non-negative.
    pub weight: f64,
}

/// Pixel value for a stellar surface-brightness map.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct StellarSurfaceBrightness {
    /// Integrated photon radiance in photons cm^-2 ns^-1 sr^-1.
    pub integrated_ph_cm2_ns_sr: f64,
    /// B-band tenth-magnitude-star surface brightness per square degree.
    pub b_s10: f64,
    /// V-band tenth-magnitude-star surface brightness per square degree.
    pub v_s10: f64,
}

impl StellarSurfaceBrightness {
    /// Return an empty pixel value.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            integrated_ph_cm2_ns_sr: 0.0,
            b_s10: 0.0,
            v_s10: 0.0,
        }
    }

    /// Return true when all values are finite and non-negative.
    #[must_use]
    pub fn is_valid(self) -> bool {
        self.integrated_ph_cm2_ns_sr.is_finite()
            && self.b_s10.is_finite()
            && self.v_s10.is_finite()
            && self.integrated_ph_cm2_ns_sr >= 0.0
            && self.b_s10 >= 0.0
            && self.v_s10 >= 0.0
    }
}

/// Provenance metadata for generated stellar maps.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StellarMapProvenance {
    /// Dataset name of the generated map.
    pub dataset_name: String,
    /// Dataset version string.
    pub version: String,
    /// UTC generation date or timestamp.
    pub generation_date_utc: String,
    /// Source catalogue name or path.
    pub source_catalogue: String,
    /// Optional source catalogue release identifier.
    pub source_catalogue_release: Option<String>,
    /// Optional source catalogue license.
    pub source_catalogue_license: Option<String>,
    /// Optional source catalogue checksum.
    pub source_catalogue_checksum: Option<String>,
    /// Optional documented magnitude limit.
    pub magnitude_limit: Option<String>,
    /// Photometric band definition.
    pub band_definition: String,
    /// Photometry model identifier.
    pub photometry_model: String,
    /// Optional smoothing description.
    pub smoothing: Option<String>,
    /// Generator identifier.
    pub generator: String,
}

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
            return Err(StellarMapError::InvalidRadianceScale(self.integrated_per_v_s10));
        }
        if self.smoothing_fwhm_deg.is_some() {
            return Err(StellarMapError::UnsupportedSmoothing);
        }

        let npix = usize::try_from(self.grid.npix()).expect("HEALPix npix fits usize");
        let mut values = vec![StellarSurfaceBrightness::zero(); npix];
        let mut included = 0_u64;

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
                values[slot].b_s10 += flux_10mag_units(mag) * record.weight;
            }
            if let Some(mag) = record.v_mag {
                values[slot].v_s10 += flux_10mag_units(mag) * record.weight;
            }
        }

        if included == 0 {
            return Err(StellarMapError::EmptyFilteredCatalogue);
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
                self.min_v_mag.map_or(true, |min| value.value() >= min.value())
                    && self.max_v_mag.map_or(true, |max| value.value() <= max.value())
            }
            None => self.min_v_mag.is_none() && self.max_v_mag.is_none(),
        }
    }
}

/// Convert an apparent magnitude to tenth-magnitude-star flux units.
#[must_use]
pub fn flux_10mag_units(magnitude: ApparentMagnitude) -> f64 {
    10.0_f64.powf(-0.4 * (magnitude.value() - 10.0))
}

/// Validate that a stellar map contains finite, non-negative values.
pub fn validate_stellar_map_values(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
) -> Result<()> {
    crate::healpix::validate_healpix_map_complete(map)?;
    if map.values.iter().all(|value| value.is_valid()) {
        Ok(())
    } else {
        Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
    }
}

/// Validate conservation of B/V flux after area-normalized binning.
pub fn validate_flux_conservation(
    input_b_flux_sum: f64,
    input_v_flux_sum: f64,
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    tolerance: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !input_b_flux_sum.is_finite() || !input_v_flux_sum.is_finite() || !tolerance.is_finite() {
        return Err(StellarMapError::Validation("flux inputs and tolerance must be finite"));
    }

    let area = map.grid.pixel_area_deg2();
    let b_sum = map.values.iter().map(|value| value.b_s10 * area).sum::<f64>();
    let v_sum = map.values.iter().map(|value| value.v_s10 * area).sum::<f64>();

    if (b_sum - input_b_flux_sum).abs() <= tolerance && (v_sum - input_v_flux_sum).abs() <= tolerance {
        Ok(())
    } else {
        Err(StellarMapError::Validation("B/V flux is not conserved"))
    }
}

/// Validate that the Galactic plane is brighter than high-latitude poles.
///
/// The validator uses a deterministic broad-brush diagnostic: mean V S10 for
/// `|b| <= 10 deg` divided by mean V S10 for `|b| >= 60 deg`.
pub fn validate_plane_pole_contrast(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    min_ratio: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !min_ratio.is_finite() || min_ratio < 0.0 {
        return Err(StellarMapError::Validation("min_ratio must be finite and non-negative"));
    }

    let mut plane_sum = 0.0;
    let mut plane_count = 0_u64;
    let mut pole_sum = 0.0;
    let mut pole_count = 0_u64;

    for (raw, value) in map.values.iter().enumerate() {
        let (_, lat) = map
            .grid
            .pixel_to_lon_lat_rad_checked(HealpixIndex::new(u64::try_from(raw).expect("index fits u64")))?;
        let lat_deg = lat.to_degrees().abs();
        if lat_deg <= 10.0 {
            plane_sum += value.v_s10;
            plane_count += 1;
        } else if lat_deg >= 60.0 {
            pole_sum += value.v_s10;
            pole_count += 1;
        }
    }

    if plane_count == 0 || pole_count == 0 {
        return Err(StellarMapError::Validation("plane/pole regions are empty"));
    }

    let plane_mean = plane_sum / plane_count as f64;
    let pole_mean = pole_sum / pole_count as f64;
    let ratio = if pole_mean == 0.0 {
        if plane_mean > 0.0 { f64::INFINITY } else { 1.0 }
    } else {
        plane_mean / pole_mean
    };

    if ratio >= min_ratio {
        Ok(())
    } else {
        Err(StellarMapError::Validation("Galactic plane/pole contrast is too low"))
    }
}

/// Validate that there is no strong artificial discontinuity at `l = 0/360`.
///
/// The diagnostic compares mean V S10 in two symmetric longitude strips near the
/// wrap seam for pixels with `|b| <= 30 deg`.
pub fn validate_no_longitude_wrap_artifact(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    max_relative_jump: f64,
) -> Result<()> {
    validate_stellar_map_values(map)?;
    if !max_relative_jump.is_finite() || max_relative_jump < 0.0 {
        return Err(StellarMapError::Validation(
            "max_relative_jump must be finite and non-negative",
        ));
    }

    let mut low_sum = 0.0;
    let mut low_count = 0_u64;
    let mut high_sum = 0.0;
    let mut high_count = 0_u64;

    for (raw, value) in map.values.iter().enumerate() {
        let (lon, lat) = map
            .grid
            .pixel_to_lon_lat_rad_checked(HealpixIndex::new(u64::try_from(raw).expect("index fits u64")))?;
        let lon_deg = lon.to_degrees();
        let lat_deg = lat.to_degrees().abs();
        if lat_deg <= 30.0 && lon_deg <= 10.0 {
            low_sum += value.v_s10;
            low_count += 1;
        } else if lat_deg <= 30.0 && lon_deg >= 350.0 {
            high_sum += value.v_s10;
            high_count += 1;
        }
    }

    if low_count == 0 || high_count == 0 {
        return Err(StellarMapError::Validation("longitude wrap regions are empty"));
    }

    let low_mean = low_sum / low_count as f64;
    let high_mean = high_sum / high_count as f64;
    let scale = low_mean.abs().max(high_mean.abs()).max(1.0);
    let jump = (low_mean - high_mean).abs() / scale;

    if jump <= max_relative_jump {
        Ok(())
    } else {
        Err(StellarMapError::Validation("longitude wrap jump exceeds tolerance"))
    }
}

/// Serialize a Galactic stellar surface-brightness map to a simple CSV string.
#[must_use]
pub fn stellar_map_to_csv(
    map: &HealpixMap<Galactic, StellarSurfaceBrightness>,
    provenance: &StellarMapProvenance,
) -> String {
    let mut out = String::new();
    out.push_str("# map_type=healpix\n");
    out.push_str(&format!("# nside={}\n", map.grid.nside().get()));
    out.push_str(&format!("# ordering={:?}\n", map.grid.ordering()).to_lowercase());
    out.push_str("# coordinate_frame=galactic\n");
    out.push_str(&format!("# dataset_name={}\n", provenance.dataset_name));
    out.push_str(&format!("# version={}\n", provenance.version));
    out.push_str(&format!("# generation_date_utc={}\n", provenance.generation_date_utc));
    out.push_str(&format!("# source_catalogue={}\n", provenance.source_catalogue));
    push_optional_header(&mut out, "source_catalogue_release", provenance.source_catalogue_release.as_deref());
    push_optional_header(&mut out, "source_catalogue_license", provenance.source_catalogue_license.as_deref());
    push_optional_header(&mut out, "source_catalogue_checksum", provenance.source_catalogue_checksum.as_deref());
    push_optional_header(&mut out, "magnitude_limit", provenance.magnitude_limit.as_deref());
    out.push_str(&format!("# photometry_model={}\n", provenance.photometry_model));
    out.push_str(&format!("# band_definition={}\n", provenance.band_definition));
    push_optional_header(&mut out, "smoothing_fwhm_deg", provenance.smoothing.as_deref());
    out.push_str("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10\n");
    for (index, value) in map.values.iter().enumerate() {
        out.push_str(&format!(
            "{},{:.17e},{:.17e},{:.17e}\n",
            index, value.integrated_ph_cm2_ns_sr, value.b_s10, value.v_s10
        ));
    }
    out
}

fn push_optional_header(out: &mut String, key: &str, value: Option<&str>) {
    if let Some(value) = value {
        out.push_str(&format!("# {key}={value}\n"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames::EquatorialMeanJ2000;
    use crate::healpix::{direction_from_lon_lat_rad, HealpixOrdering, Nside};
    use std::f64::consts::PI;

    fn provenance() -> StellarMapProvenance {
        StellarMapProvenance {
            dataset_name: "fixture".to_string(),
            version: "v1".to_string(),
            generation_date_utc: "2026-01-01T00:00:00Z".to_string(),
            source_catalogue: "deterministic fixture".to_string(),
            source_catalogue_release: Some("test".to_string()),
            source_catalogue_license: None,
            source_catalogue_checksum: None,
            magnitude_limit: Some("V<=10".to_string()),
            band_definition: "integrated 300-650 nm photon radiance plus B/V S10 diagnostics".to_string(),
            photometry_model: "v_s10_scaled_integrated_v1".to_string(),
            smoothing: None,
            generator: "siderust-test".to_string(),
        }
    }

    fn builder() -> StellarSurfaceBrightnessMapBuilder {
        StellarSurfaceBrightnessMapBuilder {
            grid: HealpixGrid::new(Nside::new(8).expect("nside"), HealpixOrdering::Ring)
                .expect("grid"),
            min_v_mag: None,
            max_v_mag: None,
            integrated_per_v_s10: 42.0,
            smoothing_fwhm_deg: None,
        }
    }

    fn galactic_record(lon_deg: f64, lat_deg: f64, v_mag: f64) -> StellarCatalogueRecord {
        let galactic: Direction<Galactic> =
            direction_from_lon_lat_rad(lon_deg.to_radians(), lat_deg.to_radians());
        let equatorial: Direction<EquatorialMeanJ2000> = galactic.to_frame();
        StellarCatalogueRecord {
            source_id: None,
            direction: equatorial,
            b_mag: Some(ApparentMagnitude::new(v_mag).expect("finite B magnitude")),
            v_mag: Some(ApparentMagnitude::new(v_mag).expect("finite V magnitude")),
            weight: 1.0,
        }
    }

    #[test]
    fn magnitude_ten_maps_to_one_tenth_mag_unit() {
        let mag = ApparentMagnitude::new(10.0).expect("finite magnitude");
        assert!((flux_10mag_units(mag) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn stellar_fixture_conserves_b_and_v_flux() {
        let records = vec![galactic_record(0.0, 0.0, 10.0), galactic_record(90.0, 0.0, 11.0)];
        let input_flux = records
            .iter()
            .map(|record| flux_10mag_units(record.v_mag.expect("V magnitude")))
            .sum::<f64>();
        let map = builder().build(records, provenance()).expect("map");
        validate_flux_conservation(input_flux, input_flux, &map, 1e-12).expect("flux conserved");
    }

    #[test]
    fn plane_concentrated_fixture_passes_plane_pole_contrast() {
        let records = vec![
            galactic_record(0.0, 0.0, 7.0),
            galactic_record(90.0, 2.0, 7.5),
            galactic_record(180.0, -3.0, 8.0),
        ];
        let map = builder().build(records, provenance()).expect("map");
        validate_plane_pole_contrast(&map, 2.0).expect("plane brighter than poles");
    }

    #[test]
    fn empty_filtered_catalogues_fail() {
        let mut cut_builder = builder();
        cut_builder.min_v_mag = Some(ApparentMagnitude::new(0.0).expect("mag"));
        cut_builder.max_v_mag = Some(ApparentMagnitude::new(1.0).expect("mag"));
        let records = vec![galactic_record(0.0, 0.0, 10.0)];
        assert!(matches!(
            cut_builder.build(records, provenance()),
            Err(StellarMapError::EmptyFilteredCatalogue)
        ));
    }

    #[test]
    fn csv_contains_expected_headers_and_rows() {
        let records = vec![galactic_record(0.0, 0.0, 10.0)];
        let provenance = provenance();
        let map = builder().build(records, provenance.clone()).expect("map");
        let csv = stellar_map_to_csv(&map, &provenance);
        assert!(csv.contains("# map_type=healpix"));
        assert!(csv.contains("# coordinate_frame=galactic"));
        assert!(csv.contains("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10"));
        assert_eq!(csv.lines().filter(|line| !line.starts_with('#')).count(), map.values.len() + 1);
    }

    #[test]
    fn longitude_wrap_validator_accepts_uniform_empty_seam() {
        let grid = HealpixGrid::new(Nside::new(8).expect("nside"), HealpixOrdering::Ring).expect("grid");
        let mut values = vec![StellarSurfaceBrightness::zero(); usize::try_from(grid.npix()).expect("npix")];
        for value in &mut values {
            value.v_s10 = 1.0;
            value.integrated_ph_cm2_ns_sr = 1.0;
        }
        let map = HealpixMap::<Galactic, StellarSurfaceBrightness>::new(grid, values).expect("map");
        validate_no_longitude_wrap_artifact(&map, 0.0).expect("no wrap jump");
    }

    #[test]
    fn full_sky_area_in_square_degrees_is_consistent() {
        let grid = builder().grid;
        let area = grid.pixel_area_deg2() * grid.npix() as f64;
        let sphere_deg2 = 4.0 * PI * (180.0 / PI).powi(2);
        assert!((area - sphere_deg2).abs() < 1e-9);
    }
}
