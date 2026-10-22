// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;
use crate::healpix::{HealpixError, HealpixIndex, HealpixMap};
use alloc::vec;

fn ring_grid(nside: u32) -> HealpixGrid {
    HealpixGrid::ring(Nside::new(nside).expect("nside")).expect("grid")
}

fn map_with_v<F>(nside: u32, brightness: F) -> HealpixMap<Galactic, StellarSurfaceBrightness>
where
    F: Fn(f64, f64) -> f64,
{
    let grid = ring_grid(nside);
    let pixels = (0..grid.npix())
        .map(|raw| {
            let center = grid
                .pixel_center_spherical::<Galactic>(HealpixIndex::new(raw))
                .expect("ring pixel center");
            StellarSurfaceBrightness {
                integrated_ph_cm2_ns_sr: 0.0,
                b_s10: 0.0,
                v_s10: brightness(center.azimuth.value(), center.polar.value()),
            }
        })
        .collect();
    HealpixMap::new(grid, pixels).expect("complete map")
}

#[test]
fn map_validation_accepts_zero_pixels_and_rejects_every_invalid_field() {
    let grid = ring_grid(1);
    let mut map = HealpixMap::<Galactic, _>::new(
        grid,
        vec![StellarSurfaceBrightness::zero(); grid.npix() as usize],
    )
    .expect("map");
    validate_stellar_map_values(&map).expect("valid zero map");

    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        map.values_mut()[0].integrated_ph_cm2_ns_sr = invalid;
        assert!(matches!(
            validate_stellar_map_values(&map),
            Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
        ));
        map.values_mut()[0] = StellarSurfaceBrightness::zero();

        map.values_mut()[0].b_s10 = invalid;
        assert!(matches!(
            validate_stellar_map_values(&map),
            Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
        ));
        map.values_mut()[0] = StellarSurfaceBrightness::zero();

        map.values_mut()[0].v_s10 = invalid;
        assert!(matches!(
            validate_stellar_map_values(&map),
            Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
        ));
        map.values_mut()[0] = StellarSurfaceBrightness::zero();
    }
    validate_stellar_map_values(&map).expect("restored map");
}

#[test]
fn flux_conservation_checks_b_and_v_independently() {
    let mut source = record(0.0, 0.0, 10.0);
    source.b_mag = Some(ApparentMagnitude::new(5.0).expect("B"));
    source.weight = 2.0;
    let map = builder().build(vec![source], provenance()).expect("map");
    let inner = map.healpix_map();

    validate_flux_conservation(200.0, 2.0, inner, 1e-10).expect("both fluxes conserved");
    assert!(matches!(
        validate_flux_conservation(201.0, 2.0, inner, 1e-10),
        Err(StellarMapError::Validation("B/V flux is not conserved"))
    ));
    assert!(matches!(
        validate_flux_conservation(200.0, 3.0, inner, 1e-10),
        Err(StellarMapError::Validation("B/V flux is not conserved"))
    ));
    validate_flux_conservation(201.0, 2.0, inner, 1.1).expect("absolute tolerance");

    let empty = map_with_v(1, |_, _| 0.0);
    validate_flux_conservation(0.0, 0.0, &empty, 0.0).expect("zero flux");
}

#[test]
fn flux_conservation_rejects_nonfinite_inputs_and_invalid_tolerances() {
    let map = map_with_v(1, |_, _| 0.0);
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            validate_flux_conservation(invalid, 0.0, &map, 0.0),
            Err(StellarMapError::Validation("flux inputs must be finite"))
        ));
        assert!(matches!(
            validate_flux_conservation(0.0, invalid, &map, 0.0),
            Err(StellarMapError::Validation("flux inputs must be finite"))
        ));
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            validate_flux_conservation(0.0, 0.0, &map, invalid),
            Err(StellarMapError::Validation("tolerance must be finite and non-negative"))
        ));
    }

    let mut bad = map_with_v(1, |_, _| 0.0);
    bad.values_mut()[0].v_s10 = f64::NAN;
    assert!(matches!(
        validate_flux_conservation(f64::NAN, 0.0, &bad, -1.0),
        Err(StellarMapError::Healpix(HealpixError::InvalidMapValue))
    ));
}

#[test]
fn plane_pole_contrast_accepts_uniform_map_at_unity_and_rejects_higher_threshold() {
    let map = map_with_v(8, |_, _| 4.0);
    validate_plane_pole_contrast(&map, 0.0).expect("zero threshold");
    validate_plane_pole_contrast(&map, 1.0).expect("uniform ratio is one");
    assert!(matches!(
        validate_plane_pole_contrast(&map, 1.1),
        Err(StellarMapError::Validation("Galactic plane/pole contrast is too low"))
    ));
}

#[test]
fn plane_pole_contrast_handles_zero_pole_and_zero_plane_flux() {
    let plane_only = map_with_v(8, |_, lat| if lat.abs() <= 10.0 { 5.0 } else { 0.0 });
    validate_plane_pole_contrast(&plane_only, f64::MAX)
        .expect("positive plane over zero pole is infinite contrast");

    let all_zero = map_with_v(8, |_, _| 0.0);
    validate_plane_pole_contrast(&all_zero, 1.0).expect("zero over zero is unity");
    assert!(matches!(
        validate_plane_pole_contrast(&all_zero, 2.0),
        Err(StellarMapError::Validation("Galactic plane/pole contrast is too low"))
    ));

    let pole_only = map_with_v(8, |_, lat| if lat.abs() >= 60.0 { 2.0 } else { 0.0 });
    validate_plane_pole_contrast(&pole_only, 0.0).expect("zero contrast meets zero threshold");
    assert!(matches!(
        validate_plane_pole_contrast(&pole_only, 0.1),
        Err(StellarMapError::Validation("Galactic plane/pole contrast is too low"))
    ));
}

#[test]
fn plane_pole_contrast_handles_empty_regions_and_invalid_thresholds() {
    let coarse_map = map_with_v(1, |_, _| 1.0);
    assert!(matches!(
        validate_plane_pole_contrast(&coarse_map, 1.0),
        Err(StellarMapError::Validation("plane/pole regions are empty"))
    ));

    let valid_map = map_with_v(8, |_, _| 1.0);
    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            validate_plane_pole_contrast(&valid_map, invalid),
            Err(StellarMapError::Validation("min_ratio must be finite and non-negative"))
        ));
    }
}

#[test]
fn longitude_seam_is_continuous_for_uniform_and_zero_maps() {
    for brightness in [0.0, 2.0] {
        let map = map_with_v(8, |_, _| brightness);
        validate_no_longitude_wrap_artifact(&map, 0.0).expect("no artificial jump");
    }
}

#[test]
fn longitude_seam_detects_asymmetric_brightness_near_zero_and_360() {
    let map = map_with_v(8, |lon, lat| {
        if lat.abs() <= 30.0 && lon <= 10.0 {
            1.0
        } else if lat.abs() <= 30.0 && lon >= 350.0 {
            9.0
        } else {
            2.0
        }
    });
    assert!(matches!(
        validate_no_longitude_wrap_artifact(&map, 0.8),
        Err(StellarMapError::Validation("longitude wrap jump exceeds tolerance"))
    ));
    validate_no_longitude_wrap_artifact(&map, 0.9)
        .expect("jump remains below relaxed threshold");
}

#[test]
fn longitude_seam_relative_jump_uses_unit_floor_for_dim_values() {
    let map = map_with_v(8, |lon, lat| {
        if lat.abs() <= 30.0 && lon <= 10.0 {
            0.2
        } else {
            0.0
        }
    });
    assert!(matches!(
        validate_no_longitude_wrap_artifact(&map, 0.15),
        Err(StellarMapError::Validation("longitude wrap jump exceeds tolerance"))
    ));
    validate_no_longitude_wrap_artifact(&map, 0.25).expect("jump is 0.2 against floor 1");
}

#[test]
fn longitude_seam_rejects_missing_strips_and_invalid_thresholds() {
    let coarse_map = map_with_v(1, |_, _| 1.0);
    assert!(matches!(
        validate_no_longitude_wrap_artifact(&coarse_map, 1.0),
        Err(StellarMapError::Validation("longitude wrap regions are empty"))
    ));

    let valid_map = map_with_v(8, |_, _| 1.0);
    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            validate_no_longitude_wrap_artifact(&valid_map, invalid),
            Err(StellarMapError::Validation("max_relative_jump must be finite and non-negative"))
        ));
    }
}

#[test]
fn spatial_validators_propagate_unsupported_healpix_ordering() {
    let grid =
        HealpixGrid::new(Nside::new(8).expect("nside"), HealpixOrdering::Nested).expect("grid");
    let map = HealpixMap::<Galactic, _>::new(
        grid,
        vec![StellarSurfaceBrightness::zero(); grid.npix() as usize],
    )
    .expect("map");
    validate_stellar_map_values(&map).expect("map values do not depend on ordering");

    assert!(matches!(
        validate_plane_pole_contrast(&map, 1.0),
        Err(StellarMapError::Healpix(HealpixError::UnsupportedOrdering(HealpixOrdering::Nested)))
    ));
    assert!(matches!(
        validate_no_longitude_wrap_artifact(&map, 1.0),
        Err(StellarMapError::Healpix(HealpixError::UnsupportedOrdering(HealpixOrdering::Nested)))
    ));
}
