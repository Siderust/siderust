use super::*;
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
fn nested_requires_power_of_two_nside() {
    let err = HealpixGrid::new(Nside::new(3).expect("valid nside"), HealpixOrdering::Nested)
        .expect_err("nested nside=3 should fail");
    assert_eq!(err, HealpixError::NestedNsideNotPowerOfTwo(3));
}
