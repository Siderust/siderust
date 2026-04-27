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

use super::{OutOfRange, TableError};

/// Validate that an axis is strictly increasing and has ≥ 2 samples.
#[inline]
pub fn validate_axis(name: &'static str, xs: &[f64]) -> Result<(), TableError> {
    if xs.len() < 2 {
        return Err(TableError::TooFewSamples { axis: name, len: xs.len() });
    }
    for i in 1..xs.len() {
        if !(xs[i] > xs[i - 1]) {
            return Err(TableError::NotMonotonic { axis: name, at_index: i });
        }
    }
    Ok(())
}

#[inline]
fn locate(xs: &[f64], x: f64) -> (usize, f64) {
    debug_assert!(xs.len() >= 2);
    // Caller should have handled out-of-range; here we assume x ∈ [xs[0], xs[-1]].
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

/// 1D linear interpolation `y(x)` over a strictly-increasing axis with the
/// given [`OutOfRange`] policy.
///
/// Mirrors [`crate::spectra::algo::interp_linear`] and is provided here so
/// table-shaped callers can stay self-contained when the `spectra` feature
/// is off.
#[inline]
pub fn linear_1d(xs: &[f64], ys: &[f64], x: f64, oor: OutOfRange) -> Result<f64, TableError> {
    debug_assert_eq!(xs.len(), ys.len());
    let (lo, hi) = (xs[0], xs[xs.len() - 1]);
    if x < lo {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(ys[0]),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error => Err(TableError::OutOfRange { axis: "x", value: x, lo, hi }),
        };
    }
    if x > hi {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(ys[ys.len() - 1]),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error => Err(TableError::OutOfRange { axis: "x", value: x, lo, hi }),
        };
    }
    let (i0, t) = locate(xs, x);
    let (y0, y1) = (ys[i0], ys[i0 + 1]);
    Ok(y0 + t * (y1 - y0))
}

/// 2D bilinear interpolation `f(x, y)` over strictly-increasing axes `xs`
/// (length `NX`) and `ys` (length `NY`), with the table laid out row-major
/// as `table[iy][ix]` (length `NY × NX`).
///
/// Per-axis [`OutOfRange`] policies are honoured independently.
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
) -> Result<f64, TableError> {
    debug_assert_eq!(xs.len(), table.first().map(|r| r.len()).unwrap_or(0));
    debug_assert_eq!(ys.len(), table.len());

    let (x_lo, x_hi) = (xs[0], xs[xs.len() - 1]);
    let (y_lo, y_hi) = (ys[0], ys[ys.len() - 1]);

    // x-axis policy
    let (x_clamped, zero_x) = if x < x_lo {
        match oor_x {
            OutOfRange::ClampToEndpoints => (x_lo, false),
            OutOfRange::Zero => (x_lo, true),
            OutOfRange::Error => {
                return Err(TableError::OutOfRange { axis: "x", value: x, lo: x_lo, hi: x_hi })
            }
        }
    } else if x > x_hi {
        match oor_x {
            OutOfRange::ClampToEndpoints => (x_hi, false),
            OutOfRange::Zero => (x_hi, true),
            OutOfRange::Error => {
                return Err(TableError::OutOfRange { axis: "x", value: x, lo: x_lo, hi: x_hi })
            }
        }
    } else {
        (x, false)
    };

    // y-axis policy
    let (y_clamped, zero_y) = if y < y_lo {
        match oor_y {
            OutOfRange::ClampToEndpoints => (y_lo, false),
            OutOfRange::Zero => (y_lo, true),
            OutOfRange::Error => {
                return Err(TableError::OutOfRange { axis: "y", value: y, lo: y_lo, hi: y_hi })
            }
        }
    } else if y > y_hi {
        match oor_y {
            OutOfRange::ClampToEndpoints => (y_hi, false),
            OutOfRange::Zero => (y_hi, true),
            OutOfRange::Error => {
                return Err(TableError::OutOfRange { axis: "y", value: y, lo: y_lo, hi: y_hi })
            }
        }
    } else {
        (y, false)
    };

    if zero_x || zero_y {
        return Ok(0.0);
    }

    let (ix0, tx) = locate(xs, x_clamped);
    let (iy0, ty) = locate(ys, y_clamped);
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
    fn linear_1d_matches_endpoint_clamp() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [10.0, 20.0, 30.0];
        assert_eq!(linear_1d(&xs, &ys, -1.0, OutOfRange::ClampToEndpoints).unwrap(), 10.0);
        assert_eq!(linear_1d(&xs, &ys, 3.0, OutOfRange::ClampToEndpoints).unwrap(), 30.0);
        assert_eq!(linear_1d(&xs, &ys, 0.5, OutOfRange::ClampToEndpoints).unwrap(), 15.0);
    }

    #[test]
    fn linear_1d_zero_policy() {
        let xs = [0.0, 1.0];
        let ys = [10.0, 20.0];
        assert_eq!(linear_1d(&xs, &ys, -1.0, OutOfRange::Zero).unwrap(), 0.0);
    }

    #[test]
    fn linear_1d_error_policy() {
        let xs = [0.0, 1.0];
        let ys = [10.0, 20.0];
        assert!(matches!(
            linear_1d(&xs, &ys, 2.0, OutOfRange::Error),
            Err(TableError::OutOfRange { .. })
        ));
    }

    #[test]
    fn bilinear_corners_recover_table_values() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [10.0, 20.0];
        let r0: &[f64] = &[1.0, 2.0, 3.0];
        let r1: &[f64] = &[10.0, 20.0, 30.0];
        let table: &[&[f64]] = &[r0, r1];
        assert_eq!(
            bilinear(&xs, &ys, table, 0.0, 10.0, OutOfRange::Error, OutOfRange::Error).unwrap(),
            1.0
        );
        assert_eq!(
            bilinear(&xs, &ys, table, 2.0, 20.0, OutOfRange::Error, OutOfRange::Error).unwrap(),
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
        let v = bilinear(&xs, &ys, table, 1.0, 1.0, OutOfRange::Error, OutOfRange::Error).unwrap();
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
        ).unwrap();
        assert_eq!(v, 0.0);
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

        let got = bilinear(&xs, &ys, t, xq, yq, OutOfRange::Error, OutOfRange::Error).unwrap();
        assert_eq!(got.to_bits(), expected.to_bits(), "bit-for-bit parity required");
    }
}
