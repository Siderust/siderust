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
//! [`Direction<F>`], not ambiguous raw `(longitude, latitude)` pairs. The raw
//! radian kernels are kept crate-private and exist only to implement the RING
//! indexing equations.
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

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::spherical;
use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::marker::PhantomData;

/// Result alias for HEALPix operations.
pub type Result<T> = std::result::Result<T, HealpixError>;

/// Error type for HEALPix grid, index, and map operations.
#[derive(Debug, thiserror::Error, Clone, PartialEq)]
pub enum HealpixError {
    /// `nside` must be greater than zero.
    #[error("HEALPix nside must be greater than zero")]
    InvalidNside,
    /// NESTED ordering requires a power-of-two `nside`.
    #[error("HEALPix NESTED ordering requires power-of-two nside, got {0}")]
    NestedNsideNotPowerOfTwo(u32),
    /// The selected ordering is not yet implemented for this operation.
    #[error("HEALPix ordering {0:?} is not supported by this operation yet")]
    UnsupportedOrdering(HealpixOrdering),
    /// A pixel index was outside the grid range.
    #[error("HEALPix pixel index {index} is outside the valid range 0..{npix}")]
    PixelIndexOutOfRange {
        /// Invalid pixel index.
        index: u64,
        /// Number of pixels in the grid.
        npix: u64,
    },
    /// The provided longitude/latitude was invalid.
    #[error("invalid HEALPix angular input: {reason}")]
    InvalidAngles {
        /// Explanation of the angular validation failure.
        reason: &'static str,
    },
    /// A map has a number of values inconsistent with its grid.
    #[error("HEALPix map length {len} does not match grid npix {npix}")]
    MapLengthMismatch {
        /// Number of stored values.
        len: usize,
        /// Number of pixels expected by the grid.
        npix: u64,
    },
    /// A finite, non-negative map value check failed.
    #[error("HEALPix map contains an invalid non-finite or negative value")]
    InvalidMapValue,
}

/// HEALPix ordering scheme.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum HealpixOrdering {
    /// HEALPix RING ordering.
    Ring,
    /// HEALPix NESTED ordering.
    Nested,
}

/// Validated HEALPix `nside` resolution parameter.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Nside(u32);

impl Nside {
    /// Construct a validated `nside`.
    pub fn new(value: u32) -> Result<Self> {
        if value == 0 {
            return Err(HealpixError::InvalidNside);
        }
        Ok(Self(value))
    }

    /// Return the underlying integer value.
    #[must_use]
    pub fn get(self) -> u32 {
        self.0
    }
}

impl TryFrom<u32> for Nside {
    type Error = HealpixError;

    fn try_from(value: u32) -> Result<Self> {
        Self::new(value)
    }
}

impl From<Nside> for u32 {
    fn from(value: Nside) -> Self {
        value.0
    }
}

/// Validated HEALPix pixel index.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct HealpixIndex(u64);

impl HealpixIndex {
    /// Construct a pixel index without grid-range validation.
    ///
    /// Use [`HealpixGrid::validate_index`] when a concrete grid is available.
    #[must_use]
    pub fn new(value: u64) -> Self {
        Self(value)
    }

    /// Return the underlying integer index.
    #[must_use]
    pub fn get(self) -> u64 {
        self.0
    }
}

impl From<HealpixIndex> for u64 {
    fn from(value: HealpixIndex) -> Self {
        value.0
    }
}

/// HEALPix grid definition.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct HealpixGrid {
    nside: Nside,
    ordering: HealpixOrdering,
}

impl HealpixGrid {
    /// Construct a validated HEALPix grid.
    pub fn new(nside: Nside, ordering: HealpixOrdering) -> Result<Self> {
        if matches!(ordering, HealpixOrdering::Nested) && !nside.get().is_power_of_two() {
            return Err(HealpixError::NestedNsideNotPowerOfTwo(nside.get()));
        }

        Ok(Self { nside, ordering })
    }

    /// Construct a RING-ordered grid.
    pub fn ring(nside: Nside) -> Result<Self> {
        Self::new(nside, HealpixOrdering::Ring)
    }

    /// Return the grid `nside`.
    #[must_use]
    pub fn nside(&self) -> Nside {
        self.nside
    }

    /// Return the grid ordering.
    #[must_use]
    pub fn ordering(&self) -> HealpixOrdering {
        self.ordering
    }

    /// Return the number of pixels, `12 * nside^2`.
    #[must_use]
    pub fn npix(&self) -> u64 {
        let nside = u64::from(self.nside.get());
        12 * nside * nside
    }

    /// Return the solid angle of each pixel in steradians.
    #[must_use]
    pub fn pixel_area_sr(&self) -> f64 {
        4.0 * PI / self.npix() as f64
    }

    /// Return the solid angle of each pixel in square degrees.
    #[must_use]
    pub fn pixel_area_deg2(&self) -> f64 {
        self.pixel_area_sr() * (180.0 / PI).powi(2)
    }

    /// Validate that `index` belongs to this grid.
    pub fn validate_index(&self, index: HealpixIndex) -> Result<()> {
        let npix = self.npix();
        if index.get() >= npix {
            return Err(HealpixError::PixelIndexOutOfRange {
                index: index.get(),
                npix,
            });
        }
        Ok(())
    }

    /// Convert a typed cartesian direction to a HEALPix pixel index.
    pub fn direction_to_pixel<F>(&self, direction: Direction<F>) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        let [x, y, z] = direction.as_array();
        let r_xy = x.hypot(y);
        let r = r_xy.hypot(z);
        if !r.is_finite() || r == 0.0 {
            return Err(HealpixError::InvalidAngles {
                reason: "direction vector must be finite and non-zero",
            });
        }
        let lon = y.atan2(x);
        let lat = (z / r).clamp(-1.0, 1.0).asin();
        self.lon_lat_rad_to_pixel_checked(lon, lat)
    }

    /// Convert a typed spherical direction to a HEALPix pixel index.
    ///
    /// Siderust spherical directions store latitude-like polar angle first and
    /// longitude-like azimuth second, both as typed angular quantities.
    pub fn spherical_to_pixel<F>(&self, direction: spherical::Direction<F>) -> Result<HealpixIndex>
    where
        F: ReferenceFrame,
    {
        self.lon_lat_rad_to_pixel_checked(
            direction.azimuth.value().to_radians(),
            direction.polar.value().to_radians(),
        )
    }

    /// Return the typed cartesian direction at the centre of a pixel.
    pub fn pixel_center<F>(&self, index: HealpixIndex) -> Result<Direction<F>>
    where
        F: ReferenceFrame,
    {
        let (lon, lat) = self.pixel_to_lon_lat_rad_checked(index)?;
        Ok(direction_from_lon_lat_rad(lon, lat))
    }

    /// Crate-private raw RING kernel from longitude/latitude in radians.
    pub(crate) fn lon_lat_rad_to_pixel_checked(
        &self,
        lon_rad: f64,
        lat_rad: f64,
    ) -> Result<HealpixIndex> {
        if self.ordering != HealpixOrdering::Ring {
            return Err(HealpixError::UnsupportedOrdering(self.ordering));
        }
        if !lon_rad.is_finite() || !lat_rad.is_finite() {
            return Err(HealpixError::InvalidAngles {
                reason: "longitude and latitude must be finite",
            });
        }
        if lat_rad < -FRAC_PI_2 || lat_rad > FRAC_PI_2 {
            return Err(HealpixError::InvalidAngles {
                reason: "latitude must lie in [-pi/2, pi/2] radians",
            });
        }

        let nside = u64::from(self.nside.get());
        let nside_f = nside as f64;
        let nl4 = 4 * nside;
        let ncap = 2 * nside * (nside - 1);
        let phi = lon_rad.rem_euclid(TAU);
        let z = lat_rad.sin();
        let za = z.abs();
        let tt = phi / FRAC_PI_2;

        let index = if za <= 2.0 / 3.0 {
            let jp = (nside_f * (0.5 + tt - 0.75 * z)).floor() as i64;
            let jm = (nside_f * (0.5 + tt + 0.75 * z)).floor() as i64;
            let ir = i64::try_from(nside).expect("nside fits i64") + 1 + jp - jm;
            let kshift = 1 - (ir & 1);
            let mut ip = (jp + jm - i64::try_from(nside).expect("nside fits i64") + kshift + 1) / 2;
            let nl4_i = i64::try_from(nl4).expect("nl4 fits i64");
            if ip > nl4_i {
                ip -= nl4_i;
            }
            if ip < 1 {
                ip += nl4_i;
            }
            ncap + u64::try_from(ir - 1).expect("equatorial ring is positive") * nl4
                + u64::try_from(ip - 1).expect("equatorial pixel is positive")
        } else {
            let tp = tt.fract();
            let tmp = nside_f * (3.0 * (1.0 - za)).sqrt();
            let jp = (tp * tmp).floor() as u64;
            let jm = ((1.0 - tp) * tmp).floor() as u64;
            let ir = jp + jm + 1;
            let mut ip = (tt * ir as f64).floor() as u64 + 1;
            let ring_pixels = 4 * ir;
            if ip > ring_pixels {
                ip -= ring_pixels;
            }
            if z >= 0.0 {
                2 * ir * (ir - 1) + ip - 1
            } else {
                self.npix() - 2 * ir * (ir + 1) + ip - 1
            }
        };

        let pixel = HealpixIndex::new(index);
        self.validate_index(pixel)?;
        Ok(pixel)
    }

    /// Crate-private raw RING kernel to longitude/latitude in radians.
    pub(crate) fn pixel_to_lon_lat_rad_checked(
        &self,
        index: HealpixIndex,
    ) -> Result<(f64, f64)> {
        if self.ordering != HealpixOrdering::Ring {
            return Err(HealpixError::UnsupportedOrdering(self.ordering));
        }
        self.validate_index(index)?;

        let pix = index.get();
        let nside = u64::from(self.nside.get());
        let nside_f = nside as f64;
        let nl4 = 4 * nside;
        let ncap = 2 * nside * (nside - 1);
        let equatorial_end = ncap + nl4 * (2 * nside + 1);
        let fact1 = 1.0 / (1.5 * nside_f);
        let fact2 = 1.0 / (3.0 * nside_f * nside_f);

        let (phi, z) = if pix < ncap {
            let iring = ((1.0 + (1.0 + 2.0 * pix as f64).sqrt()) * 0.5).floor() as u64;
            let iphi = pix + 1 - 2 * iring * (iring - 1);
            let z = 1.0 - iring as f64 * iring as f64 * fact2;
            let phi = (iphi as f64 - 0.5) * FRAC_PI_2 / iring as f64;
            (phi, z)
        } else if pix < equatorial_end {
            let ip = pix - ncap;
            let iring = ip / nl4 + nside;
            let iphi = ip % nl4 + 1;
            let fodd = if ((iring + nside) & 1) == 1 { 1.0 } else { 0.5 };
            let z = (2.0 * nside_f - iring as f64) * fact1;
            let phi = (iphi as f64 - fodd) * FRAC_PI_2 / nside_f;
            (phi, z)
        } else {
            let ip = self.npix() - pix;
            let iring = ((1.0 + (2.0 * ip as f64 - 1.0).sqrt()) * 0.5).floor() as u64;
            let iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
            let z = -1.0 + iring as f64 * iring as f64 * fact2;
            let phi = (iphi as f64 - 0.5) * FRAC_PI_2 / iring as f64;
            (phi, z)
        };

        Ok((phi.rem_euclid(TAU), z.clamp(-1.0, 1.0).asin()))
    }
}

/// Frame-typed HEALPix map container.
#[derive(Debug, Clone, PartialEq)]
pub struct HealpixMap<F, T>
where
    F: ReferenceFrame,
{
    /// HEALPix grid shared by all values.
    pub grid: HealpixGrid,
    /// Pixel values in the grid ordering.
    pub values: Vec<T>,
    frame: PhantomData<F>,
}

impl<F, T> HealpixMap<F, T>
where
    F: ReferenceFrame,
{
    /// Construct a map and validate that `values.len() == grid.npix()`.
    pub fn new(grid: HealpixGrid, values: Vec<T>) -> Result<Self> {
        let len = values.len();
        let npix = grid.npix();
        if u64::try_from(len).expect("usize length fits u64") != npix {
            return Err(HealpixError::MapLengthMismatch { len, npix });
        }
        Ok(Self {
            grid,
            values,
            frame: PhantomData,
        })
    }

    /// Return the frame marker carried by this map type.
    #[must_use]
    pub fn frame_marker(&self) -> PhantomData<F> {
        self.frame
    }
}

/// Validate that a HEALPix map has one value for every pixel.
pub fn validate_healpix_map_complete<F, T>(map: &HealpixMap<F, T>) -> Result<()>
where
    F: ReferenceFrame,
{
    let len = map.values.len();
    let npix = map.grid.npix();
    if u64::try_from(len).expect("usize length fits u64") != npix {
        return Err(HealpixError::MapLengthMismatch { len, npix });
    }
    Ok(())
}

/// Build a typed cartesian direction from longitude/latitude in radians.
pub(crate) fn direction_from_lon_lat_rad<F>(lon_rad: f64, lat_rad: f64) -> Direction<F>
where
    F: ReferenceFrame,
{
    let cos_lat = lat_rad.cos();
    Direction::<F>::from_array([cos_lat * lon_rad.cos(), cos_lat * lon_rad.sin(), lat_rad.sin()])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames;
    use crate::qtty::Degrees;

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
            let direction: Direction<frames::Galactic> = grid.pixel_center(index).expect("pixel centre");
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
    fn spherical_direction_api_is_frame_typed() {
        let grid = ring_grid(4);
        let direction = spherical::Direction::<frames::Galactic>::new_unchecked(
            Degrees::new(0.0),
            Degrees::new(0.0),
        );
        let pixel = grid.spherical_to_pixel(direction).expect("pixel");
        grid.validate_index(pixel).expect("valid pixel");
    }

    #[test]
    fn nested_requires_power_of_two_nside() {
        let err = HealpixGrid::new(Nside::new(3).expect("valid nside"), HealpixOrdering::Nested)
            .expect_err("nested nside=3 should fail");
        assert_eq!(err, HealpixError::NestedNsideNotPowerOfTwo(3));
    }
}
