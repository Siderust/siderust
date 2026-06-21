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
//! transforms, bins them into a typed
//! [`crate::healpix::HealpixMap<crate::coordinates::frames::Galactic, StellarSurfaceBrightness>`],
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

mod brightness;
mod builder;
mod csv;
mod error;
mod magnitude;
mod photometry;
mod provenance;
mod record;
mod validation;

pub use brightness::StellarSurfaceBrightness;
pub use builder::StellarSurfaceBrightnessMapBuilder;
pub use csv::stellar_map_to_csv;
pub use error::{Result, StellarMapError};
pub use magnitude::ApparentMagnitude;
pub use photometry::flux_10mag_units;
pub use provenance::StellarMapProvenance;
pub use record::StellarCatalogueRecord;
pub use validation::{
    validate_flux_conservation, validate_no_longitude_wrap_artifact,
    validate_plane_pole_contrast, validate_stellar_map_values,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::cartesian::Direction;
    use crate::coordinates::frames::{EquatorialMeanJ2000, Galactic};
    use crate::coordinates::transform::TransformFrame;
    use crate::healpix::{
        direction_from_lon_lat_rad, HealpixGrid, HealpixMap, HealpixOrdering, Nside,
    };
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
            band_definition: "integrated 300-650 nm photon radiance plus B/V S10 diagnostics"
                .to_string(),
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
        let records = vec![
            galactic_record(0.0, 0.0, 10.0),
            galactic_record(90.0, 0.0, 11.0),
        ];
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
    fn catalogues_without_usable_photometry_fail() {
        let mut record = galactic_record(0.0, 0.0, 10.0);
        record.b_mag = None;
        record.v_mag = None;
        assert!(matches!(
            builder().build(vec![record], provenance()),
            Err(StellarMapError::NoUsablePhotometry)
        ));
    }

    #[test]
    fn csv_contains_expected_headers_and_rows() {
        let records = vec![galactic_record(0.0, 0.0, 10.0)];
        let provenance = provenance();
        let map = builder().build(records, provenance.clone()).expect("map");
        let csv = stellar_map_to_csv(&map, &provenance);
        assert!(csv.contains("# map_type=healpix"));
        assert!(csv.contains("# schema_version=stellar_healpix_csv_v1"));
        assert!(csv.contains("# coordinate_frame=galactic"));
        assert!(csv.contains("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10"));
        assert_eq!(
            csv.lines().filter(|line| !line.starts_with('#')).count(),
            map.values().len() + 1
        );
    }

    #[test]
    fn longitude_wrap_validator_accepts_uniform_empty_seam() {
        let grid = HealpixGrid::new(Nside::new(8).expect("nside"), HealpixOrdering::Ring)
            .expect("grid");
        let mut values =
            vec![StellarSurfaceBrightness::zero(); usize::try_from(grid.npix()).expect("npix")];
        for value in &mut values {
            value.v_s10 = 1.0;
            value.integrated_ph_cm2_ns_sr = 1.0;
        }
        let map =
            HealpixMap::<Galactic, StellarSurfaceBrightness>::new(grid, values).expect("map");
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
