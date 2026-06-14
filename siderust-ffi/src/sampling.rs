// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Sky-grid sampling exposed over the C ABI.
//!
//! Mirrors [`siderust::coordinates::SkyGrid`] / [`SkyGridCell`]: a hemispherical
//! alt/az grid sampler.  A single flattening entry point materialises every
//! cell as `(altitude, azimuth, solid_angle)` into a heap-allocated C array
//! that the caller frees with [`siderust_sky_grid_cells_free`].

use siderust::coordinates::{SkyGrid, SkyGridCell};
use siderust::qtty::Degrees;

use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c};

/// A single sky-grid cell: a Horizontal direction and its approximate solid
/// angle.
///
/// Mirrors [`siderust::coordinates::SkyGridCell`].
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiderustSkyGridCell {
    /// Altitude above the horizon (polar component), in degrees, `[-90, 90]`.
    pub altitude_deg: f64,
    /// Azimuth from North (clockwise), in degrees, `[0, 360)`.
    pub azimuth_deg: f64,
    /// Approximate solid angle subtended by the cell, in steradians.
    pub solid_angle_sr: f64,
}

impl SiderustSkyGridCell {
    #[inline]
    fn ffi_from(cell: &SkyGridCell) -> Self {
        Self {
            altitude_deg: cell.direction.polar.value(),
            azimuth_deg: cell.direction.azimuth.value(),
            solid_angle_sr: cell.solid_angle.value(),
        }
    }
}

/// Materialise every cell of a hemispherical alt/az grid.
///
/// The grid covers the altitude band `[alt_min_deg, alt_max_deg)` with altitude
/// step `alt_step_deg` and azimuth step `az_step_deg`.  When `equal_area` is
/// true, the azimuth count per altitude ring scales with `cos(alt)` (the
/// [`SkyGrid::equal_area`] construction); otherwise a fixed azimuth step is used
/// ([`SkyGrid::with_steps`]).
///
/// # Safety
/// `out` and `count` must be valid, writable pointers.  On success `*out` owns
/// a heap array of `*count` cells that must be released with
/// [`siderust_sky_grid_cells_free`].
#[no_mangle]
pub unsafe extern "C" fn siderust_sky_grid_cells(
    alt_min_deg: f64,
    alt_max_deg: f64,
    alt_step_deg: f64,
    az_step_deg: f64,
    equal_area: bool,
    out: *mut *mut SiderustSkyGridCell,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() || count.is_null() {
            return SiderustStatus::NullPointer;
        }
        let base = if equal_area {
            SkyGrid::equal_area(Degrees::new(alt_step_deg), Degrees::new(az_step_deg))
        } else {
            SkyGrid::with_steps(Degrees::new(alt_step_deg), Degrees::new(az_step_deg))
        };
        let grid = base.with_alt_range(Degrees::new(alt_min_deg), Degrees::new(alt_max_deg));
        let cells: Vec<SkyGridCell> = grid.iter_cells().collect();
        vec_to_c(cells, SiderustSkyGridCell::ffi_from, out, count)
    }}
}

/// Free a sky-grid cell array produced by [`siderust_sky_grid_cells`].
///
/// # Safety
/// `ptr` and `count` must originate from the same `siderust_sky_grid_cells`
/// call and must not have been freed before; `ptr` must not be used afterwards.
#[no_mangle]
pub unsafe extern "C" fn siderust_sky_grid_cells_free(ptr: *mut SiderustSkyGridCell, count: usize) {
    unsafe { free_boxed_slice(ptr, count) };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    #[test]
    fn uniform_grid_has_expected_cell_count() {
        let mut ptr_out: *mut SiderustSkyGridCell = ptr::null_mut();
        let mut count: usize = 0;
        // Uniform 10° grid over [0,90): 9 altitude rings × 36 azimuth steps.
        let status = unsafe {
            siderust_sky_grid_cells(0.0, 90.0, 10.0, 10.0, false, &mut ptr_out, &mut count)
        };
        assert_eq!(status, SiderustStatus::Ok);
        assert_eq!(count, 9 * 36);
        assert!(!ptr_out.is_null());
        unsafe { siderust_sky_grid_cells_free(ptr_out, count) };
    }

    #[test]
    fn equal_area_solid_angles_sum_to_hemisphere() {
        let mut ptr_out: *mut SiderustSkyGridCell = ptr::null_mut();
        let mut count: usize = 0;
        let status =
            unsafe { siderust_sky_grid_cells(0.0, 90.0, 5.0, 5.0, true, &mut ptr_out, &mut count) };
        assert_eq!(status, SiderustStatus::Ok);
        let cells = unsafe { std::slice::from_raw_parts(ptr_out, count) };
        let total: f64 = cells.iter().map(|c| c.solid_angle_sr).sum();
        let hemisphere = 2.0 * std::f64::consts::PI;
        assert!((total - hemisphere).abs() / hemisphere < 0.01);
        unsafe { siderust_sky_grid_cells_free(ptr_out, count) };
    }

    #[test]
    fn null_out_is_rejected() {
        let mut count: usize = 0;
        let status = unsafe {
            siderust_sky_grid_cells(0.0, 90.0, 10.0, 10.0, false, ptr::null_mut(), &mut count)
        };
        assert_eq!(status, SiderustStatus::NullPointer);
    }
}
