// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;
use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames;
use crate::coordinates::spherical;
use std::f64::consts::PI;

fn ring_grid(nside: u32) -> HealpixGrid {
    HealpixGrid::ring(Nside::new(nside).expect("valid nside")).expect("valid grid")
}

#[test]
fn nside_pixel_counts_match_definition() {
    assert_eq!(ring_grid(1).npix(), 12);
    assert_eq!(ring_grid(2).npix(), 48);
    assert_eq!(ring_grid(64).npix(), 49_152);
}

#[test]
fn nside_rejects_overflowing_pixel_counts() {
    assert_eq!(
        Nside::new(MAX_NSIDE + 1),
        Err(HealpixError::NsideTooLarge(MAX_NSIDE + 1))
    );
}

#[test]
fn pixel_areas_sum_to_full_sphere() {
    let grid = ring_grid(64);
    let area = grid.pixel_area_sr() * grid.npix() as f64;
    assert!((area - 4.0 * PI).abs() < 1e-12);
}

#[test]
fn pixel_centres_roundtrip_through_direction_api() {
    let grid = ring_grid(8);
    for raw in [0, 1, 7, 63, 128, 511, grid.npix() - 1] {
        let index = HealpixIndex::new(raw);
        let direction: Direction<frames::Galactic> = grid.pixel_center(index).expect("center");
        assert_eq!(grid.direction_to_pixel(direction).expect("pixel"), index);
    }
}

#[test]
fn spherical_pixel_centres_match_cartesian_centres_and_roundtrip() {
    let grid = ring_grid(8);
    for raw in [0, 1, 7, 63, 128, 511, grid.npix() - 1] {
        let index = HealpixIndex::new(raw);
        let cartesian: Direction<frames::Galactic> = grid.pixel_center(index).expect("center");
        let spherical: spherical::Direction<frames::Galactic> = grid
            .pixel_center_spherical(index)
            .expect("spherical center");

        let from_spherical = spherical.to_cartesian();
        for (actual, expected) in from_spherical.as_array().iter().zip(cartesian.as_array()) {
            assert!((actual - expected).abs() < 1e-12);
        }
        assert_eq!(
            grid.direction_to_pixel(from_spherical).expect("pixel"),
            index
        );
    }
}

#[test]
fn spherical_pixel_centres_are_canonical_at_wrap_and_polar_regions() {
    let grid = ring_grid(8);
    for raw in [0, 1, grid.npix() / 2, grid.npix() - 2, grid.npix() - 1] {
        let center: spherical::Direction<frames::Galactic> = grid
            .pixel_center_spherical(HealpixIndex::new(raw))
            .expect("spherical center");
        assert!((-90.0..=90.0).contains(&center.polar.value()));
        assert!((0.0..360.0).contains(&center.azimuth.value()));
    }
}

#[test]
fn spherical_pixel_centres_preserve_nested_ordering_behavior() {
    let grid = HealpixGrid::new(Nside::new(8).expect("valid nside"), HealpixOrdering::Nested)
        .expect("valid grid");
    assert!(matches!(
        grid.pixel_center_spherical::<frames::Galactic>(HealpixIndex::new(0)),
        Err(HealpixError::UnsupportedOrdering(HealpixOrdering::Nested))
    ));
}

#[test]
fn nested_requires_power_of_two_nside() {
    let err = HealpixGrid::new(Nside::new(3).expect("valid nside"), HealpixOrdering::Nested)
        .expect_err("nested nside=3 should fail");
    assert_eq!(err, HealpixError::NestedNsideNotPowerOfTwo(3));
}
