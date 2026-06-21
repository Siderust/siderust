use super::*;
use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{EquatorialMeanJ2000, Galactic};
use crate::coordinates::transform::TransformFrame;
use crate::healpix::{direction_from_theta_phi, HealpixGrid, HealpixOrdering, Nside};
use std::f64::consts::FRAC_PI_2;

fn provenance() -> StellarMapProvenance {
    StellarMapProvenance {
        dataset_name: "fixture".to_string(),
        version: "v1".to_string(),
        generation_date_utc: "2026-01-01T00:00:00Z".to_string(),
        source_catalogue: "fixture".to_string(),
        source_catalogue_release: None,
        source_catalogue_license: None,
        source_catalogue_checksum: None,
        magnitude_limit: Some("V<=10".to_string()),
        band_definition: "B/V S10".to_string(),
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
    }
}

fn record(phi_deg: f64, dec_deg: f64, v_mag: f64) -> StellarCatalogueRecord {
    let theta = FRAC_PI_2 - dec_deg.to_radians();
    let phi = phi_deg.to_radians();
    let galactic: Direction<Galactic> = direction_from_theta_phi(theta, phi);
    let equatorial: Direction<EquatorialMeanJ2000> = galactic.to_frame();
    StellarCatalogueRecord {
        source_id: None,
        direction: equatorial,
        b_mag: Some(ApparentMagnitude::new(v_mag).expect("B")),
        v_mag: Some(ApparentMagnitude::new(v_mag).expect("V")),
        weight: 1.0,
    }
}

#[test]
fn magnitude_ten_maps_to_one_flux_unit() {
    let mag = ApparentMagnitude::new(10.0).expect("magnitude");
    assert!((flux_10mag_units(mag) - 1.0).abs() < 1e-15);
}

#[test]
fn generated_map_keeps_provenance() {
    let provenance = provenance();
    let map = builder()
        .build(vec![record(0.0, 0.0, 10.0)], provenance.clone())
        .expect("map");
    assert_eq!(map.provenance(), &provenance);
}

#[test]
fn catalogues_without_usable_photometry_fail() {
    let mut source = record(0.0, 0.0, 10.0);
    source.b_mag = None;
    source.v_mag = None;
    assert!(matches!(
        builder().build(vec![source], provenance()),
        Err(StellarMapError::NoUsablePhotometry)
    ));
}

#[test]
fn csv_contains_header_and_one_row_per_pixel() {
    let map = builder()
        .build(vec![record(0.0, 0.0, 10.0)], provenance())
        .expect("map");
    let csv = csv::to_csv(&map);
    assert!(csv.starts_with("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10"));
    assert_eq!(csv.lines().count(), map.values().len() + 1);
}
