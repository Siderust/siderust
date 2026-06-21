// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! HEALPix full-sky grid primitives.
//!
//! ## Scientific scope
//!
//! HEALPix partitions the celestial sphere into equal-area pixels. It is widely
//! used for full-sky maps because each pixel has the same solid angle and the
//! RING ordering has efficient latitude-band structure.
//!
//! ## Technical scope
//!
//! This module implements the typed Siderust-facing API for HEALPix grids and
//! maps. Public callers work with frame-typed astronomical directions such as
//! [`crate::coordinates::cartesian::Direction`], not ambiguous raw
//! `(longitude, latitude)` pairs. The raw radian kernels are internal
//! implementation details in the private `ring` module.
//!
//! The first supported ordering is [`HealpixOrdering::Ring`].
//! [`HealpixOrdering::Nested`] is validated for power-of-two `nside` values but
//! pixel conversions currently return [`HealpixError::UnsupportedOrdering`].
//!
//! ## References
//!
//! - Gorski, K. M. et al. (2005). *HEALPix: A Framework for High-Resolution
//!   Discretization and Fast Analysis of Data Distributed on the Sphere*, ApJ,
//!   622, 759.

mod error;
mod grid;
mod index;
mod map;
mod nside;
mod ordering;
mod ring;
mod validation;

pub use error::{HealpixError, Result};
pub use grid::HealpixGrid;
pub use index::HealpixIndex;
pub use map::HealpixMap;
pub use nside::{Nside, MAX_NSIDE};
pub use ordering::HealpixOrdering;
pub use validation::validate_healpix_map_complete;

#[cfg(test)]
pub(crate) use ring::direction_from_lon_lat_rad;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::cartesian::Direction;
    use crate::coordinates::frames;
    use crate::coordinates::spherical;
    use crate::qtty::Degrees;
    use std::f64::consts::{PI, TAU};

    fn ring_grid(nside: u32) -> HealpixGrid {
        HealpixGrid::ring(Nside::new(nside).expect("valid nside")).expect("valid ring grid")
    }

    #[test]
    fn nside_pixel_counts_match_healpix_definition() {
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
    fn pixel_centres_are_finite_and_roundtrip_for_representative_indices() {
        let grid = ring_grid(8);
        let indices = [0, 1, 7, 63, 128, 511, grid.npix() - 1];

        for raw in indices {
            let index = HealpixIndex::new(raw);
            let direction: Direction<frames::Galactic> =
                grid.pixel_center(index).expect("pixel centre");
            let [x, y, z] = direction.as_array();
            assert!(x.is_finite());
            assert!(y.is_finite());
            assert!(z.is_finite());
            assert_eq!(grid.direction_to_pixel(direction).expect("roundtrip"), index);
        }
    }

    #[test]
    fn longitude_wrap_inputs_are_valid_and_continuous() {
        let grid = ring_grid(16);
        let a = grid
            .lon_lat_rad_to_pixel_checked((-1.0e-12_f64).rem_euclid(TAU), 0.0)
            .expect("valid wrapped lon");
        let b = grid
            .lon_lat_rad_to_pixel_checked(TAU - 1.0e-12_f64, 0.0)
            .expect("valid high lon");
        assert_eq!(a, b);
    }

    #[test]
    fn spherical_direction_api_matches_cartesian_for_nonzero_latitude() {
        let grid = ring_grid(16);
        let spherical = spherical::Direction::<frames::Galactic>::new_unchecked(
            Degrees::new(30.0),
            Degrees::new(45.0),
        );
        let cartesian: Direction<frames::Galactic> =
            direction_from_lon_lat_rad(45.0_f64.to_radians(), 30.0_f64.to_radians());

        assert_eq!(
            grid.spherical_to_pixel(spherical).expect("spherical pixel"),
            grid.direction_to_pixel(cartesian).expect("cartesian pixel")
        );
    }

    #[test]
    fn nested_requires_power_of_two_nside() {
        let err = HealpixGrid::new(Nside::new(3).expect("valid nside"), HealpixOrdering::Nested)
            .expect_err("nested nside=3 should fail");
        assert_eq!(err, HealpixError::NestedNsideNotPowerOfTwo(3));
    }
}
