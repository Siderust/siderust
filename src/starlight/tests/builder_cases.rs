// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;
use crate::healpix::{HealpixError, HealpixMap};
use alloc::format;
use alloc::vec;
use alloc::vec::Vec;

fn magnitude(value: f64) -> ApparentMagnitude {
    ApparentMagnitude::new(value).expect("finite magnitude")
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "expected {expected}, got {actual}"
    );
}

#[test]
fn builder_rejects_invalid_radiance_scales() {
    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut settings = builder();
        settings.integrated_per_v_s10 = invalid;
        assert!(matches!(
            settings.build(vec![record(0.0, 0.0, 10.0)], provenance()),
            Err(StellarMapError::InvalidRadianceScale(value))
                if (value.is_nan() && invalid.is_nan()) || value == invalid
        ));
    }
}

#[test]
fn builder_rejects_invalid_weights_even_for_excluded_sources() {
    let mut settings = builder();
    settings.min_v_mag = Some(magnitude(20.0));

    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut source = record(0.0, 0.0, 10.0);
        source.weight = invalid;
        assert!(matches!(
            settings.build(vec![source], provenance()),
            Err(StellarMapError::InvalidWeight(value))
                if (value.is_nan() && invalid.is_nan()) || value == invalid
        ));
    }
}

#[test]
fn builder_distinguishes_empty_input_filtered_input_and_missing_photometry() {
    assert!(matches!(
        builder().build(vec![], provenance()),
        Err(StellarMapError::EmptyFilteredCatalogue)
    ));

    let mut settings = builder();
    settings.min_v_mag = Some(magnitude(11.0));
    assert!(matches!(
        settings.build(vec![record(0.0, 0.0, 10.0)], provenance()),
        Err(StellarMapError::EmptyFilteredCatalogue)
    ));

    let mut missing = record(0.0, 0.0, 10.0);
    missing.b_mag = None;
    missing.v_mag = None;
    assert!(matches!(
        builder().build(vec![missing], provenance()),
        Err(StellarMapError::NoUsablePhotometry)
    ));

    let mut zero_weight = record(0.0, 0.0, 10.0);
    zero_weight.weight = 0.0;
    assert!(matches!(
        builder().build(vec![zero_weight], provenance()),
        Err(StellarMapError::NoUsablePhotometry)
    ));

    let tiny_flux = record(0.0, 0.0, 1.0e6);
    assert!(matches!(
        builder().build(vec![tiny_flux], provenance()),
        Err(StellarMapError::NoUsablePhotometry)
    ));
}

#[test]
fn magnitude_bounds_are_inclusive_and_evaluated_on_v_only() {
    let mut settings = builder();
    settings.min_v_mag = Some(magnitude(9.0));
    settings.max_v_mag = Some(magnitude(11.0));

    let mut b_only = record(0.0, 0.0, 10.0);
    b_only.v_mag = None;
    b_only.b_mag = Some(magnitude(-10.0));

    let map = settings
        .build(
            vec![
                record(0.0, 0.0, 8.0),
                record(0.0, 0.0, 9.0),
                record(0.0, 0.0, 10.0),
                record(0.0, 0.0, 11.0),
                record(0.0, 0.0, 12.0),
                b_only,
            ],
            provenance(),
        )
        .expect("inclusive bounds");

    let area = map.grid().pixel_area_deg2();
    let expected_flux = [9.0, 10.0, 11.0]
        .into_iter()
        .map(|v| flux_10mag_units(magnitude(v)))
        .sum::<f64>();
    let actual_flux = map.values().iter().map(|p| p.v_s10 * area).sum::<f64>();
    assert_close(actual_flux, expected_flux);

    settings.min_v_mag = Some(magnitude(12.0));
    settings.max_v_mag = Some(magnitude(9.0));
    assert!(matches!(
        settings.build(vec![record(0.0, 0.0, 10.0)], provenance()),
        Err(StellarMapError::EmptyFilteredCatalogue)
    ));
}

#[test]
fn single_sided_magnitude_bounds_filter_independently() {
    let mut min_only = builder();
    min_only.min_v_mag = Some(magnitude(10.0));
    let selected = min_only
        .build(
            vec![record(0.0, 0.0, 9.0), record(0.0, 0.0, 10.0)],
            provenance(),
        )
        .expect("lower bound includes endpoint");
    let area = selected.grid().pixel_area_deg2();
    assert_close(selected.values().iter().map(|p| p.v_s10 * area).sum(), 1.0);

    let mut max_only = builder();
    max_only.max_v_mag = Some(magnitude(10.0));
    let selected = max_only
        .build(
            vec![record(0.0, 0.0, 10.0), record(0.0, 0.0, 11.0)],
            provenance(),
        )
        .expect("upper bound includes endpoint");
    let area = selected.grid().pixel_area_deg2();
    assert_close(selected.values().iter().map(|p| p.v_s10 * area).sum(), 1.0);
}

#[test]
fn b_only_source_is_usable_without_a_v_filter_but_not_with_one() {
    let mut source = record(0.0, 0.0, 10.0);
    source.v_mag = None;

    let map = builder()
        .build(vec![source.clone()], provenance())
        .expect("B-only photometry");
    let area = map.grid().pixel_area_deg2();
    assert_close(map.values().iter().map(|p| p.b_s10 * area).sum(), 1.0);
    assert!(map
        .values()
        .iter()
        .all(|p| p.v_s10 == 0.0 && p.integrated_ph_cm2_ns_sr == 0.0));

    for (min, max) in [(Some(magnitude(0.0)), None), (None, Some(magnitude(20.0)))] {
        let mut settings = builder();
        settings.min_v_mag = min;
        settings.max_v_mag = max;
        assert!(matches!(
            settings.build(vec![source.clone()], provenance()),
            Err(StellarMapError::EmptyFilteredCatalogue)
        ));
    }
}

#[test]
fn builder_accumulates_b_and_v_independently_with_weights_and_area() {
    let mut first = record(42.0, -20.0, 10.0);
    first.b_mag = Some(magnitude(9.0));
    first.weight = 2.0;

    let mut second = record(42.0, -20.0, 15.0);
    second.b_mag = None;
    second.weight = 3.0;

    let map = builder()
        .build(vec![first, second], provenance())
        .expect("weighted map");

    let populated: Vec<_> = map.values().iter().filter(|p| p.v_s10 > 0.0).collect();
    assert_eq!(populated.len(), 1, "coincident stars must share a pixel");

    let area = map.grid().pixel_area_deg2();
    let pixel = populated[0];
    assert_close(pixel.b_s10 * area, 2.0 * flux_10mag_units(magnitude(9.0)));
    assert_close(
        pixel.v_s10 * area,
        2.0 + 3.0 * flux_10mag_units(magnitude(15.0)),
    );
    assert_close(pixel.integrated_ph_cm2_ns_sr, 42.0 * pixel.v_s10);
    assert!(map.values().iter().all(|p| p.is_valid()));
}

#[test]
fn v_only_and_b_only_sources_contribute_to_their_respective_bands() {
    let mut v_only = record(0.0, 0.0, 15.0);
    v_only.b_mag = None;
    let mut b_only = record(0.0, 0.0, 10.0);
    b_only.v_mag = None;

    let map = builder()
        .build(vec![v_only, b_only], provenance())
        .expect("partial bands");
    let area = map.grid().pixel_area_deg2();
    assert_close(map.values().iter().map(|p| p.b_s10 * area).sum(), 1.0);
    assert_close(map.values().iter().map(|p| p.v_s10 * area).sum(), 0.01);
}

#[test]
fn zero_radiance_calibration_keeps_s10_values() {
    let mut settings = builder();
    settings.integrated_per_v_s10 = 0.0;
    let map = settings
        .build(vec![record(0.0, 0.0, 10.0)], provenance())
        .expect("zero calibration");
    assert!(map.values().iter().any(|p| p.v_s10 > 0.0));
    assert!(map
        .values()
        .iter()
        .all(|p| p.integrated_ph_cm2_ns_sr == 0.0));
}

#[test]
fn builder_rejects_overflowed_photometry_or_integrated_radiance() {
    let mut bright = record(0.0, 0.0, -1_000.0);
    bright.b_mag = None;
    assert!(matches!(
        builder().build(vec![bright], provenance()),
        Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
    ));

    let mut settings = builder();
    settings.integrated_per_v_s10 = f64::MAX;
    assert!(matches!(
        settings.build(vec![record(0.0, 0.0, -100.0)], provenance()),
        Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
    ));
}

#[test]
fn nested_ordering_is_reported_as_unsupported_by_binning() {
    let mut settings = builder();
    settings.grid = HealpixGrid::new(Nside::new(8).expect("nside"), HealpixOrdering::Nested)
        .expect("nested grid");
    assert!(matches!(
        settings.build(vec![record(0.0, 0.0, 10.0)], provenance()),
        Err(StellarMapError::Healpix(HealpixError::UnsupportedOrdering(
            HealpixOrdering::Nested
        )))
    ));
}

#[test]
fn map_accessors_and_into_parts_preserve_grid_pixels_and_provenance() {
    let metadata = provenance();
    let map = builder()
        .build(vec![record(0.0, 0.0, 10.0)], metadata.clone())
        .expect("map");
    assert_eq!(map.grid(), builder().grid);
    assert_eq!(map.healpix_map().grid(), map.grid());
    assert_eq!(map.healpix_map().values(), map.values());

    let original_values = map.values().to_vec();
    let (inner, extracted_metadata) = map.into_parts();
    assert_eq!(inner.values(), original_values);
    assert_eq!(inner.grid(), builder().grid);
    assert_eq!(extracted_metadata, metadata);
}

#[test]
fn csv_rows_follow_healpix_order_and_preserve_numeric_values() {
    let mut settings = builder();
    settings.grid = HealpixGrid::ring(Nside::new(1).expect("nside")).expect("ring");
    let map = settings
        .build(vec![record(0.0, 0.0, 10.0)], provenance())
        .expect("map");
    let csv = csv::to_csv(&map);
    let mut lines = csv.lines();
    assert_eq!(
        lines.next(),
        Some("healpix_index,integrated_ph_cm2_ns_sr,b_s10,v_s10")
    );
    for (index, (line, pixel)) in lines.zip(map.values()).enumerate() {
        assert_eq!(
            line,
            format!(
                "{index},{},{},{}",
                pixel.integrated_ph_cm2_ns_sr, pixel.b_s10, pixel.v_s10
            )
        );
    }
    assert_eq!(csv.lines().count(), 13);
    assert_eq!(csv.matches(",0,0,0\n").count(), 11);
}

#[test]
fn manually_created_map_can_be_wrapped_with_provenance() {
    let grid = HealpixGrid::ring(Nside::new(1).expect("nside")).expect("grid");
    let pixels = vec![StellarSurfaceBrightness::zero(); grid.npix() as usize];
    let inner = HealpixMap::<Galactic, _>::new(grid, pixels.clone()).expect("map");
    let metadata = provenance();
    let wrapped = StellarSurfaceBrightnessMap::new(inner, metadata.clone());
    assert_eq!(wrapped.values(), pixels);
    assert_eq!(wrapped.provenance(), &metadata);
}
