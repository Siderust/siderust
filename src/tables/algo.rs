// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Untyped `f64` interpolation kernels for [`crate::tables`].
//!
//! These primitives operate on raw slices and `&[&[f64]]`-like views and
//! exist for callers that need bit-for-bit `numpy.interp`-style parity with
//! existing pipelines (NSB's Leinert lookup, ESO Sky-Model tables, …).
//!
//! Typed callers should prefer [`Grid1D::interp_at`](super::Grid1D::interp_at)
//! and [`Grid2D::interp_at`](super::Grid2D::interp_at).

use super::{AxisDirection, OutOfRange, TableError};

/// Validate that an axis is strictly monotonic (either direction) and has ≥ 2
/// samples. Returns the detected [`AxisDirection`].
///
/// Equal consecutive values are rejected with [`TableError::NotMonotonic`].
#[inline]
pub fn validate_axis(name: &'static str, xs: &[f64]) -> Result<AxisDirection, TableError> {
    if xs.len() < 2 {
        return Err(TableError::TooFewSamples { axis: name, len: xs.len() });
    }
    // Determine direction from the first step, then verify consistency.
    let dir = if xs[1] > xs[0] {
        AxisDirection::Ascending
    } else if xs[1] < xs[0] {
        AxisDirection::Descending
    } else {
        return Err(TableError::NotMonotonic { axis: name, at_index: 1 });
    };
    for i in 2..xs.len() {
        let ok = match dir {
            AxisDirection::Ascending => xs[i] > xs[i - 1],
            AxisDirection::Descending => xs[i] < xs[i - 1],
        };
        if !ok {
            return Err(TableError::NotMonotonic { axis: name, at_index: i });
        }
    }
    Ok(dir)
}

/// Returns `(lo, hi)` as *value-space* range, i.e. `lo ≤ hi` regardless of
/// axis direction.
#[inline]
fn axis_range(xs: &[f64], dir: AxisDirection) -> (f64, f64) {
    match dir {
        AxisDirection::Ascending => (xs[0], xs[xs.len() - 1]),
        AxisDirection::Descending => (xs[xs.len() - 1], xs[0]),
    }
}

/// Apply an out-of-range policy for a single axis.
///
/// `lo ≤ hi` must hold (use [`axis_range`] to obtain them). Returns
/// `(clamped_value, is_zero)` or an error.
#[inline]
fn apply_oor(
    x: f64,
    lo: f64,
    hi: f64,
    oor: OutOfRange,
    axis: &'static str,
) -> Result<(f64, bool), TableError> {
    if x < lo {
        match oor {
            OutOfRange::ClampToEndpoints => Ok((lo, false)),
            OutOfRange::Zero => Ok((lo, true)),
            OutOfRange::Error => {
                Err(TableError::OutOfRange { axis, value: x, lo, hi })
            }
        }
    } else if x > hi {
        match oor {
            OutOfRange::ClampToEndpoints => Ok((hi, false)),
            OutOfRange::Zero => Ok((hi, true)),
            OutOfRange::Error => {
                Err(TableError::OutOfRange { axis, value: x, lo, hi })
            }
        }
    } else {
        Ok((x, false))
    }
}

#[inline]
fn locate(xs: &[f64], x: f64, dir: AxisDirection) -> (usize, f64) {
    debug_assert!(xs.len() >= 2);
    // Caller must have clamped x to the axis range; here we assume
    // x ∈ [min(xs), max(xs)].
    match dir {
        AxisDirection::Ascending => {
            if x <= xs[0] {
                return (0, 0.0);
            }
            if x >= xs[xs.len() - 1] {
                return (xs.len() - 2, 1.0);
            }
            let i = xs.partition_point(|&v| v <= x);
            let i0 = i - 1;
            let t = (x - xs[i0]) / (xs[i0 + 1] - xs[i0]);
            (i0, t)
        }
        AxisDirection::Descending => {
            // xs[0] is the largest value, xs[n-1] is the smallest.
            if x >= xs[0] {
                return (0, 0.0);
            }
            if x <= xs[xs.len() - 1] {
                return (xs.len() - 2, 1.0);
            }
            // First index where xs[i] ≤ x (i.e. xs goes below x).
            let i = xs.partition_point(|&v| v > x);
            let i0 = i - 1;
            // Both numerator and denominator are negative; t ∈ (0, 1).
            let t = (x - xs[i0]) / (xs[i0 + 1] - xs[i0]);
            (i0, t)
        }
    }
}

/// 1D linear interpolation `y(x)` over a strictly-monotonic axis with the
/// given [`OutOfRange`] policy.
///
/// Accepts axes in either direction (pass the [`AxisDirection`] returned by
/// [`validate_axis`]). For a descending axis `ClampToEndpoints` clamps to
/// `xs[0]` / `ys[0]` when `x > xs[0]` and to `xs[n-1]` / `ys[n-1]` when
/// `x < xs[n-1]`.
///
/// Mirrors [`crate::spectra::algo::interp_linear`] and is provided here so
/// table-shaped callers can stay self-contained when the `spectra` feature
/// is off.
#[inline]
pub fn linear_1d(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
    dir: AxisDirection,
) -> Result<f64, TableError> {
    debug_assert_eq!(xs.len(), ys.len());
    let (lo, hi) = axis_range(xs, dir);
    let (lo_y, hi_y) = match dir {
        AxisDirection::Ascending => (ys[0], ys[ys.len() - 1]),
        AxisDirection::Descending => (ys[ys.len() - 1], ys[0]),
    };
    if x < lo {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(lo_y),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error => Err(TableError::OutOfRange { axis: "x", value: x, lo, hi }),
        };
    }
    if x > hi {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(hi_y),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error => Err(TableError::OutOfRange { axis: "x", value: x, lo, hi }),
        };
    }
    let (i0, t) = locate(xs, x, dir);
    let (y0, y1) = (ys[i0], ys[i0 + 1]);
    Ok(y0 + t * (y1 - y0))
}

/// 2D bilinear interpolation `f(x, y)` over strictly-monotonic axes `xs`
/// (length `NX`) and `ys` (length `NY`), with the table laid out row-major
/// as `table[iy][ix]` (length `NY × NX`).
///
/// Axes may be ascending or descending; pass the [`AxisDirection`] returned
/// by [`validate_axis`] for each. Per-axis [`OutOfRange`] policies are
/// honoured independently.
///
/// **Numerical convention** (matters for bit-for-bit parity): the interpolant
/// is computed as
///
/// ```text
/// row_lo = T[iy0][ix0] + tx · (T[iy0][ix1] - T[iy0][ix0])
/// row_hi = T[iy1][ix0] + tx · (T[iy1][ix1] - T[iy1][ix0])
/// out    = row_lo      + ty · (row_hi      - row_lo)
/// ```
///
/// — i.e. interpolate along `x` first within each row, then along `y`. This
/// matches the ordering used by NSB's `leinert_lookup_s10` so callers
/// migrating from inline math get bit-identical results.
#[inline]
pub fn bilinear(
    xs: &[f64],
    ys: &[f64],
    table: &[&[f64]],
    x: f64,
    y: f64,
    oor_x: OutOfRange,
    oor_y: OutOfRange,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
) -> Result<f64, TableError> {
    debug_assert_eq!(xs.len(), table.first().map(|r| r.len()).unwrap_or(0));
    debug_assert_eq!(ys.len(), table.len());

    let (x_lo, x_hi) = axis_range(xs, dir_x);
    let (y_lo, y_hi) = axis_range(ys, dir_y);

    let (x_clamped, zero_x) = apply_oor(x, x_lo, x_hi, oor_x, "x")?;
    let (y_clamped, zero_y) = apply_oor(y, y_lo, y_hi, oor_y, "y")?;

    if zero_x || zero_y {
        return Ok(0.0);
    }

    let (ix0, tx) = locate(xs, x_clamped, dir_x);
    let (iy0, ty) = locate(ys, y_clamped, dir_y);
    let ix1 = ix0 + 1;
    let iy1 = iy0 + 1;

    let row_lo = table[iy0][ix0] + tx * (table[iy0][ix1] - table[iy0][ix0]);
    let row_hi = table[iy1][ix0] + tx * (table[iy1][ix1] - table[iy1][ix0]);
    Ok(row_lo + ty * (row_hi - row_lo))
}

/// Bilinear interpolation kernel from four already-located corner values
/// and two pre-computed unit fractions, using the same operation ordering
/// as [`bilinear`].
///
/// Useful when the caller has its own non-uniform / wrapped indexing math
/// (e.g. the Leinert zodiacal lookup) but still wants bit-for-bit parity
/// with the shared kernel.
///
/// Corners are named after standard table indexing: `f00 = T[iy0][ix0]`,
/// `f10 = T[iy0][ix1]`, `f01 = T[iy1][ix0]`, `f11 = T[iy1][ix1]`.
/// `tx` and `ty` are unit fractions, typically in `[0, 1]`.
#[inline]
pub fn bilinear_unit(f00: f64, f10: f64, f01: f64, f11: f64, tx: f64, ty: f64) -> f64 {
    let row_lo = f00 + tx * (f10 - f00);
    let row_hi = f01 + tx * (f11 - f01);
    row_lo + ty * (row_hi - row_lo)
}

/// 3D trilinear interpolation `f(x, y, z)` over strictly-monotonic axes
/// `xs` (length `NX`), `ys` (length `NY`), `zs` (length `NZ`), with the
/// table laid out as a flat slice of length `NX·NY·NZ` using the storage
/// convention `index = (iz·NY + iy)·NX + ix` (innermost stride is `x`,
/// then `y`, then `z`).
///
/// Axes may be ascending or descending; pass the [`AxisDirection`] returned
/// by [`validate_axis`] for each. Per-axis [`OutOfRange`] policies are
/// honoured independently.
///
/// **Numerical convention** (matters for bit-for-bit parity): the interpolant
/// is computed as
///
/// ```text
/// plane_lo = bilinear_unit(f000, f100, f010, f110, tx, ty)
/// plane_hi = bilinear_unit(f001, f101, f011, f111, tx, ty)
/// out      = plane_lo + tz · (plane_hi − plane_lo)
/// ```
///
/// — i.e. interpolate along `x` first, then `y`, then `z`.
#[inline]
pub fn trilinear(
    xs: &[f64],
    ys: &[f64],
    zs: &[f64],
    table: &[f64],
    nx: usize,
    ny: usize,
    x: f64,
    y: f64,
    z: f64,
    oor_x: OutOfRange,
    oor_y: OutOfRange,
    oor_z: OutOfRange,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
    dir_z: AxisDirection,
) -> Result<f64, TableError> {
    debug_assert_eq!(table.len(), nx * ny * zs.len());

    let (x_lo, x_hi) = axis_range(xs, dir_x);
    let (y_lo, y_hi) = axis_range(ys, dir_y);
    let (z_lo, z_hi) = axis_range(zs, dir_z);

    let (x_clamped, zero_x) = apply_oor(x, x_lo, x_hi, oor_x, "x")?;
    let (y_clamped, zero_y) = apply_oor(y, y_lo, y_hi, oor_y, "y")?;
    let (z_clamped, zero_z) = apply_oor(z, z_lo, z_hi, oor_z, "z")?;

    if zero_x || zero_y || zero_z {
        return Ok(0.0);
    }

    let (ix0, tx) = locate(xs, x_clamped, dir_x);
    let (iy0, ty) = locate(ys, y_clamped, dir_y);
    let (iz0, tz) = locate(zs, z_clamped, dir_z);
    let ix1 = ix0 + 1;
    let iy1 = iy0 + 1;
    let iz1 = iz0 + 1;

    // index = (iz * NY + iy) * NX + ix
    let idx = |iz: usize, iy: usize, ix: usize| (iz * ny + iy) * nx + ix;
    let f000 = table[idx(iz0, iy0, ix0)];
    let f100 = table[idx(iz0, iy0, ix1)];
    let f010 = table[idx(iz0, iy1, ix0)];
    let f110 = table[idx(iz0, iy1, ix1)];
    let f001 = table[idx(iz1, iy0, ix0)];
    let f101 = table[idx(iz1, iy0, ix1)];
    let f011 = table[idx(iz1, iy1, ix0)];
    let f111 = table[idx(iz1, iy1, ix1)];

    Ok(trilinear_unit(f000, f100, f010, f110, f001, f101, f011, f111, tx, ty, tz))
}

/// Trilinear interpolation kernel from eight already-located corner values
/// and three pre-computed unit fractions, mirroring [`bilinear_unit`].
///
/// Corner naming follows `f_xyz` with bit 0 = x, bit 1 = y, bit 2 = z:
/// `f000 = T[iz0][iy0][ix0]`, `f111 = T[iz1][iy1][ix1]`, etc.
/// `tx`, `ty`, `tz` are unit fractions, typically in `[0, 1]`.
///
/// Implementation order (same as [`trilinear`]):
///
/// ```text
/// plane_lo = bilinear_unit(f000, f100, f010, f110, tx, ty)
/// plane_hi = bilinear_unit(f001, f101, f011, f111, tx, ty)
/// out      = plane_lo + tz · (plane_hi − plane_lo)
/// ```
#[inline]
pub fn trilinear_unit(
    f000: f64,
    f100: f64,
    f010: f64,
    f110: f64,
    f001: f64,
    f101: f64,
    f011: f64,
    f111: f64,
    tx: f64,
    ty: f64,
    tz: f64,
) -> f64 {
    let plane_lo = bilinear_unit(f000, f100, f010, f110, tx, ty);
    let plane_hi = bilinear_unit(f001, f101, f011, f111, tx, ty);
    plane_lo + tz * (plane_hi - plane_lo)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_rejects_short_axis() {
        assert!(matches!(
            validate_axis("x", &[1.0]),
            Err(TableError::TooFewSamples { .. })
        ));
    }

    #[test]
    fn validate_rejects_non_monotonic() {
        assert!(matches!(
            validate_axis("x", &[1.0, 2.0, 2.0]),
            Err(TableError::NotMonotonic { at_index: 2, .. })
        ));
    }

    #[test]
    fn validate_accepts_descending() {
        assert_eq!(
            validate_axis("x", &[3.0, 2.0, 1.0]).unwrap(),
            AxisDirection::Descending
        );
    }

    #[test]
    fn validate_accepts_ascending() {
        assert_eq!(
            validate_axis("x", &[1.0, 2.0, 3.0]).unwrap(),
            AxisDirection::Ascending
        );
    }

    #[test]
    fn validate_rejects_non_monotonic_after_descent() {
        // Goes down then up — not strictly monotonic.
        assert!(matches!(
            validate_axis("x", &[3.0, 2.0, 2.5]),
            Err(TableError::NotMonotonic { at_index: 2, .. })
        ));
    }

    #[test]
    fn linear_1d_matches_endpoint_clamp() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [10.0, 20.0, 30.0];
        assert_eq!(
            linear_1d(&xs, &ys, -1.0, OutOfRange::ClampToEndpoints, AxisDirection::Ascending)
                .unwrap(),
            10.0
        );
        assert_eq!(
            linear_1d(&xs, &ys, 3.0, OutOfRange::ClampToEndpoints, AxisDirection::Ascending)
                .unwrap(),
            30.0
        );
        assert_eq!(
            linear_1d(&xs, &ys, 0.5, OutOfRange::ClampToEndpoints, AxisDirection::Ascending)
                .unwrap(),
            15.0
        );
    }

    #[test]
    fn linear_1d_zero_policy() {
        let xs = [0.0, 1.0];
        let ys = [10.0, 20.0];
        assert_eq!(
            linear_1d(&xs, &ys, -1.0, OutOfRange::Zero, AxisDirection::Ascending).unwrap(),
            0.0
        );
    }

    #[test]
    fn linear_1d_error_policy() {
        let xs = [0.0, 1.0];
        let ys = [10.0, 20.0];
        assert!(matches!(
            linear_1d(&xs, &ys, 2.0, OutOfRange::Error, AxisDirection::Ascending),
            Err(TableError::OutOfRange { .. })
        ));
    }

    /// Descending axis interpolation gives identical results to the
    /// axis-and-values-reversed ascending equivalent.
    #[test]
    fn linear_1d_descending_matches_ascending_equivalent() {
        let xs_desc = [3.0, 2.0, 1.0];
        let ys_desc = [30.0, 20.0, 10.0];

        // Query x=2.5: between xs_desc[0]=3 and xs_desc[1]=2.
        let v = linear_1d(&xs_desc, &ys_desc, 2.5, OutOfRange::ClampToEndpoints,
                          AxisDirection::Descending).unwrap();
        assert_eq!(v, 25.0);

        // Same query on the ascending flip must give bit-identical results.
        let xs_asc = [1.0, 2.0, 3.0];
        let ys_asc = [10.0, 20.0, 30.0];
        let v_asc = linear_1d(&xs_asc, &ys_asc, 2.5, OutOfRange::ClampToEndpoints,
                              AxisDirection::Ascending).unwrap();
        assert_eq!(v.to_bits(), v_asc.to_bits(), "bit-for-bit parity required");
    }

    #[test]
    fn linear_1d_descending_clamp_low_end() {
        // x below the range [1,3] → clamp to xs[n-1]=1 → ys[n-1]=10
        let xs = [3.0, 2.0, 1.0];
        let ys = [30.0, 20.0, 10.0];
        assert_eq!(
            linear_1d(&xs, &ys, 0.0, OutOfRange::ClampToEndpoints, AxisDirection::Descending)
                .unwrap(),
            10.0
        );
    }

    #[test]
    fn linear_1d_descending_clamp_high_end() {
        // x above the range [1,3] → clamp to xs[0]=3 → ys[0]=30
        let xs = [3.0, 2.0, 1.0];
        let ys = [30.0, 20.0, 10.0];
        assert_eq!(
            linear_1d(&xs, &ys, 4.0, OutOfRange::ClampToEndpoints, AxisDirection::Descending)
                .unwrap(),
            30.0
        );
    }

    #[test]
    fn bilinear_corners_recover_table_values() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [10.0, 20.0];
        let r0: &[f64] = &[1.0, 2.0, 3.0];
        let r1: &[f64] = &[10.0, 20.0, 30.0];
        let table: &[&[f64]] = &[r0, r1];
        assert_eq!(
            bilinear(&xs, &ys, table, 0.0, 10.0, OutOfRange::Error, OutOfRange::Error,
                     AxisDirection::Ascending, AxisDirection::Ascending).unwrap(),
            1.0
        );
        assert_eq!(
            bilinear(&xs, &ys, table, 2.0, 20.0, OutOfRange::Error, OutOfRange::Error,
                     AxisDirection::Ascending, AxisDirection::Ascending).unwrap(),
            30.0
        );
    }

    #[test]
    fn bilinear_midpoint_is_average_of_four_corners() {
        let xs = [0.0, 2.0];
        let ys = [0.0, 2.0];
        let r0: &[f64] = &[1.0, 3.0];
        let r1: &[f64] = &[5.0, 7.0];
        let table: &[&[f64]] = &[r0, r1];
        let v = bilinear(&xs, &ys, table, 1.0, 1.0, OutOfRange::Error, OutOfRange::Error,
                         AxisDirection::Ascending, AxisDirection::Ascending).unwrap();
        assert_eq!(v, (1.0 + 3.0 + 5.0 + 7.0) / 4.0);
    }

    #[test]
    fn bilinear_clamp_outside_returns_corner() {
        let xs = [0.0, 1.0];
        let ys = [0.0, 1.0];
        let r0: &[f64] = &[1.0, 2.0];
        let r1: &[f64] = &[3.0, 4.0];
        let table: &[&[f64]] = &[r0, r1];
        let v = bilinear(
            &xs, &ys, table, -1.0, -1.0,
            OutOfRange::ClampToEndpoints, OutOfRange::ClampToEndpoints,
            AxisDirection::Ascending, AxisDirection::Ascending,
        ).unwrap();
        assert_eq!(v, 1.0);
    }

    #[test]
    fn bilinear_zero_policy_outside_yields_zero() {
        let xs = [0.0, 1.0];
        let ys = [0.0, 1.0];
        let r0: &[f64] = &[1.0, 2.0];
        let r1: &[f64] = &[3.0, 4.0];
        let table: &[&[f64]] = &[r0, r1];
        let v = bilinear(
            &xs, &ys, table, -1.0, 0.5,
            OutOfRange::Zero, OutOfRange::ClampToEndpoints,
            AxisDirection::Ascending, AxisDirection::Ascending,
        ).unwrap();
        assert_eq!(v, 0.0);
    }

    /// `bilinear` on a table with a descending y-axis produces bit-for-bit
    /// identical results compared to `bilinear_unit` with manually computed
    /// descending-axis indices (replicating NSB's leinert_lookup_s10 logic).
    #[test]
    fn bilinear_descending_y_matches_bilinear_unit_bit_for_bit() {
        // 4-column (x/β ascending) × 4-row (y/λ descending) synthetic table.
        // xs = β axis:  [0, 5, 10, 15] (ascending)
        // ys = λ axis: [15, 10, 5, 0] (descending, like Leinert row layout)
        // table[iy][ix] = (iy+1)*100 + (ix+1)*10
        let xs: [f64; 4] = [0.0, 5.0, 10.0, 15.0];
        let ys: [f64; 4] = [15.0, 10.0, 5.0, 0.0];
        let r0: [f64; 4] = [110.0, 120.0, 130.0, 140.0]; // ys[0]=15
        let r1: [f64; 4] = [210.0, 220.0, 230.0, 240.0]; // ys[1]=10
        let r2: [f64; 4] = [310.0, 320.0, 330.0, 340.0]; // ys[2]=5
        let r3: [f64; 4] = [410.0, 420.0, 430.0, 440.0]; // ys[3]=0
        let table: &[&[f64]] = &[&r0, &r1, &r2, &r3];

        // Query: x=7.5, y=7.5
        let xq = 7.5_f64;
        let yq = 7.5_f64;

        // Manual descending-y indexing (mirrors NSB leinert_lookup_s10):
        // ix0 = floor(xq/5) = 1, bt = (7.5 - 5)/5 = 0.5
        // ys[1]=10 >= yq=7.5 > ys[2]=5 → iy0 = 1
        // lt = (yq - ys[1]) / (ys[2] - ys[1]) = (7.5 - 10) / (5 - 10) = 0.5
        let ix0 = 1usize;
        let iy0 = 1usize;
        let bt = (xq - xs[ix0]) / (xs[ix0 + 1] - xs[ix0]);
        let lt = (yq - ys[iy0]) / (ys[iy0 + 1] - ys[iy0]);
        let expected = bilinear_unit(
            r1[ix0], r1[ix0 + 1],
            r2[ix0], r2[ix0 + 1],
            bt, lt,
        );

        let got = bilinear(
            &xs, &ys, table, xq, yq,
            OutOfRange::Error, OutOfRange::Error,
            AxisDirection::Ascending, AxisDirection::Descending,
        ).unwrap();

        assert_eq!(got.to_bits(), expected.to_bits(),
                   "bilinear(descending y) must be bit-for-bit equal to bilinear_unit: \
                    got={got}, expected={expected}");
    }

    #[test]
    fn trilinear_unit_midpoint() {
        // f(x,y,z) = 100x + 10y + z on the unit cube [0,1]^3
        // corners at (x,y,z) in {0,1}^3:
        let f000 = 0.0_f64;   // (0,0,0)
        let f100 = 100.0;     // (1,0,0)
        let f010 = 10.0;      // (0,1,0)
        let f110 = 110.0;     // (1,1,0)
        let f001 = 1.0;       // (0,0,1)
        let f101 = 101.0;     // (1,0,1)
        let f011 = 11.0;      // (0,1,1)
        let f111 = 111.0;     // (1,1,1)
        // At midpoint (0.5, 0.5, 0.5): expected = 100*0.5 + 10*0.5 + 0.5 = 55.5
        let v = trilinear_unit(f000, f100, f010, f110, f001, f101, f011, f111,
                               0.5, 0.5, 0.5);
        assert!((v - 55.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn trilinear_corners_recover_values() {
        // 2×2×2 grid; f(x,y,z) = 100x + 10y + z
        let xs = [0.0_f64, 1.0];
        let ys = [0.0_f64, 1.0];
        let zs = [0.0_f64, 1.0];
        // table[(iz*NY+iy)*NX+ix]
        let table: Vec<f64> = (0..8)
            .map(|i| {
                let ix = i & 1;
                let iy = (i >> 1) & 1;
                let iz = (i >> 2) & 1;
                (100 * ix + 10 * iy + iz) as f64
            })
            .collect();
        for (xv, yv, zv) in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)] {
            let expected = 100.0 * xv + 10.0 * yv + zv;
            let got = trilinear(&xs, &ys, &zs, &table, 2, 2,
                                xv, yv, zv,
                                OutOfRange::Error, OutOfRange::Error, OutOfRange::Error,
                                AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
                .unwrap();
            assert!((got - expected).abs() < 1e-12, "at ({xv},{yv},{zv}): got {got}, expected {expected}");
        }
    }

    #[test]
    fn trilinear_midpoint_333_grid() {
        // 3×3×3 grid; f(x,y,z) = 100x + 10y + z
        let xs = [0.0_f64, 1.0, 2.0];
        let ys = [0.0_f64, 1.0, 2.0];
        let zs = [0.0_f64, 1.0, 2.0];
        let mut table = vec![0.0_f64; 27];
        for iz in 0..3 {
            for iy in 0..3 {
                for ix in 0..3 {
                    table[(iz * 3 + iy) * 3 + ix] =
                        100.0 * ix as f64 + 10.0 * iy as f64 + iz as f64;
                }
            }
        }
        let v = trilinear(&xs, &ys, &zs, &table, 3, 3,
                          0.5, 0.5, 0.5,
                          OutOfRange::Error, OutOfRange::Error, OutOfRange::Error,
                          AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
            .unwrap();
        assert!((v - 55.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn trilinear_oor_policies() {
        let xs = [0.0_f64, 1.0];
        let ys = [0.0_f64, 1.0];
        let zs = [0.0_f64, 1.0];
        let table: Vec<f64> = (0..8)
            .map(|i| {
                let ix = i & 1;
                let iy = (i >> 1) & 1;
                let iz = (i >> 2) & 1;
                (100 * ix + 10 * iy + iz) as f64
            })
            .collect();
        // Clamp x below → same as x=0
        let v = trilinear(&xs, &ys, &zs, &table, 2, 2,
                          -1.0, 0.5, 0.5,
                          OutOfRange::ClampToEndpoints, OutOfRange::ClampToEndpoints, OutOfRange::ClampToEndpoints,
                          AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
            .unwrap();
        let expected = trilinear(&xs, &ys, &zs, &table, 2, 2,
                                 0.0, 0.5, 0.5,
                                 OutOfRange::Error, OutOfRange::Error, OutOfRange::Error,
                                 AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
            .unwrap();
        assert_eq!(v.to_bits(), expected.to_bits());

        // Zero policy for z-axis OOR
        let v = trilinear(&xs, &ys, &zs, &table, 2, 2,
                          0.5, 0.5, 5.0,
                          OutOfRange::ClampToEndpoints, OutOfRange::ClampToEndpoints, OutOfRange::Zero,
                          AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
            .unwrap();
        assert_eq!(v, 0.0);

        // Error policy
        let r = trilinear(&xs, &ys, &zs, &table, 2, 2,
                          0.5, 5.0, 0.5,
                          OutOfRange::Error, OutOfRange::Error, OutOfRange::Error,
                          AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending);
        assert!(matches!(r, Err(TableError::OutOfRange { axis: "y", .. })));
    }

    #[test]
    fn trilinear_unit_matches_trilinear_on_2x2x2() {
        // Verify trilinear_unit and trilinear give bit-for-bit identical results.
        let xs = [0.0_f64, 1.0];
        let ys = [0.0_f64, 1.0];
        let zs = [0.0_f64, 1.0];
        let table: Vec<f64> = (0..8)
            .map(|i| {
                let ix = i & 1;
                let iy = (i >> 1) & 1;
                let iz = (i >> 2) & 1;
                (100 * ix + 10 * iy + iz + 1) as f64
            })
            .collect();
        let (xq, yq, zq) = (0.3, 0.6, 0.8);
        let via_trilinear = trilinear(&xs, &ys, &zs, &table, 2, 2,
                                     xq, yq, zq,
                                     OutOfRange::Error, OutOfRange::Error, OutOfRange::Error,
                                     AxisDirection::Ascending, AxisDirection::Ascending, AxisDirection::Ascending)
            .unwrap();
        let idx = |iz: usize, iy: usize, ix: usize| (iz * 2 + iy) * 2 + ix;
        let via_unit = trilinear_unit(
            table[idx(0,0,0)], table[idx(0,0,1)],
            table[idx(0,1,0)], table[idx(0,1,1)],
            table[idx(1,0,0)], table[idx(1,0,1)],
            table[idx(1,1,0)], table[idx(1,1,1)],
            xq, yq, zq,
        );
        assert_eq!(via_trilinear.to_bits(), via_unit.to_bits(),
                   "trilinear_unit must match trilinear bit-for-bit");
    }

    /// Bit-for-bit parity check against NSB's hand-rolled bilinear in
    /// `nsb::components::zodiacal::leinert_lookup_s10`. Replays its exact
    /// numerical sequence on a small synthetic table.
    #[test]
    fn bilinear_matches_nsb_ordering_bit_for_bit() {
        // Synthetic 3×4 table; arbitrary values exercising all four corners.
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [0.0, 1.0, 2.0];
        let r0: &[f64] = &[10.5, 11.5, 12.5, 13.5];
        let r1: &[f64] = &[20.5, 21.5, 22.5, 23.5];
        let r2: &[f64] = &[30.5, 31.5, 32.5, 33.5];
        let t: &[&[f64]] = &[r0, r1, r2];

        let xq = 1.3_f64;
        let yq = 0.7_f64;

        // Compute reference using exact NSB ordering:
        //   bx0 = T[iy0][ix0] + tx*(T[iy0][ix1] - T[iy0][ix0])
        //   bx1 = T[iy1][ix0] + tx*(T[iy1][ix1] - T[iy1][ix0])
        //   out = bx0 + ty*(bx1 - bx0)
        let ix0 = 1usize;
        let iy0 = 0usize;
        let tx = (xq - xs[ix0]) / (xs[ix0 + 1] - xs[ix0]);
        let ty = (yq - ys[iy0]) / (ys[iy0 + 1] - ys[iy0]);
        let bx0 = t[iy0][ix0] + tx * (t[iy0][ix0 + 1] - t[iy0][ix0]);
        let bx1 = t[iy0 + 1][ix0] + tx * (t[iy0 + 1][ix0 + 1] - t[iy0 + 1][ix0]);
        let expected = bx0 + ty * (bx1 - bx0);

        let got = bilinear(&xs, &ys, t, xq, yq, OutOfRange::Error, OutOfRange::Error,
                           AxisDirection::Ascending, AxisDirection::Ascending).unwrap();
        assert_eq!(got.to_bits(), expected.to_bits(), "bit-for-bit parity required");
    }
}
