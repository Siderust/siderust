// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed 2D bilinearly-interpolated table.

use crate::ext_qtty::{Quantity, Scalar, Unit};
use crate::interp::OutOfRange;
use crate::provenance::Provenance;

use super::{algo, AxisDirection, TableError};

/// Half-open rectangular region in the **query coordinate system** that, if
/// the query `(x, y)` falls inside it, short-circuits the bilinear lookup
/// and returns a constant value.
///
/// Bounds are interpreted as:
///   - `x_min_inclusive`: `x >= x_min_inclusive` (when `Some`); always satisfied otherwise.
///   - `x_max_exclusive`: `x <  x_max_exclusive` (when `Some`); always satisfied otherwise.
///   - `y_min_inclusive` / `y_max_exclusive`: same semantics for `y`.
///
/// Use [`ConstantRegion::lower_corner`] for the common
/// `if x < x_max && y < y_max` pattern that motivated this API
/// (NSB Leinert zodiacal corner clamps).
#[derive(Debug, Clone, Copy)]
pub struct ConstantRegion<S: Scalar> {
    pub x_min_inclusive: Option<S>,
    pub x_max_exclusive: Option<S>,
    pub y_min_inclusive: Option<S>,
    pub y_max_exclusive: Option<S>,
    pub value: S,
}

impl<S: Scalar> ConstantRegion<S> {
    /// Constant fill region for `x < x_max && y < y_max`.
    pub fn lower_corner(x_max_exclusive: S, y_max_exclusive: S, value: S) -> Self {
        Self {
            x_min_inclusive: None,
            x_max_exclusive: Some(x_max_exclusive),
            y_min_inclusive: None,
            y_max_exclusive: Some(y_max_exclusive),
            value,
        }
    }

    #[inline]
    fn contains(&self, x: S, y: S) -> bool {
        if let Some(b) = self.x_min_inclusive { if !(x >= b) { return false; } }
        if let Some(b) = self.x_max_exclusive { if !(x <  b) { return false; } }
        if let Some(b) = self.y_min_inclusive { if !(y >= b) { return false; } }
        if let Some(b) = self.y_max_exclusive { if !(y <  b) { return false; } }
        true
    }
}

/// Strictly-monotonic 2D table of values laid out row-major as
/// `table[iy][ix]`. The two axes carry their own typed [`Unit`]
/// markers (`X` for the column axis, `Y` for the row axis); the cell
/// values carry a third [`Unit`] marker `V`.
///
/// `S` defaults to `f64`. The interpolant matches the ordering used by
/// [`algo::bilinear`] (interpolate along `x` first within each row, then
/// along `y`).
///
/// The grid optionally carries:
///   - a list of [`ConstantRegion`]s checked **before** interpolation in the
///     query-coordinate system, and
///   - a *y-axis reflection*: callers may construct the grid with
///     [`Grid2D::from_raw_row_major_y_descending`], pass the table in its
///     natural row-descending form, and still query in the same descending
///     coordinate system; internally the storage is normalized to ascending
///     and the query `y` is reflected as `y_internal = (y_max + y_min) − y`.
///     This guarantees bit-for-bit equivalence with the legacy
///     `(y_max − y − step·i0) / step` arithmetic that NSB's hand-rolled
///     Leinert lookup used.
#[derive(Debug, Clone)]
pub struct Grid2D<X: Unit, Y: Unit, V: Unit, S: Scalar = f64> {
    xs: Vec<S>,
    ys: Vec<S>,
    /// Row-major storage: row `iy` starts at `iy * NX`.
    table: Vec<S>,
    nx: usize,
    ny: usize,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
    /// When `Some(off)`, queries are interpreted in a descending coordinate
    /// system; the y-value passed to the kernel is `off − y_user`. `off`
    /// equals `ys_user[0] + ys_user[ny-1]` of the original (descending)
    /// axis supplied at construction time.
    y_reflect_offset: Option<S>,
    regions: Vec<ConstantRegion<S>>,
    provenance: Provenance,
    _markers: core::marker::PhantomData<(X, Y, V)>,
}

impl<X: Unit, Y: Unit, V: Unit, S: Scalar + Into<f64> + From<f64>>
    Grid2D<X, Y, V, S>
{
    /// Build a `Grid2D` from raw scalar slices already expressed in the
    /// declared `X` / `Y` / `V` units. The two axes are validated to be
    /// strictly monotonic; the table is required to be `ys.len() · xs.len()`
    /// long, laid out row-major.
    pub fn from_raw_row_major(
        xs: Vec<S>,
        ys: Vec<S>,
        table: Vec<S>,
    ) -> Result<Self, TableError> {
        let nx = xs.len();
        let ny = ys.len();
        if table.len() != nx * ny {
            return Err(TableError::ShapeMismatch {
                expected_x: nx,
                expected_y: ny,
                actual_rows: if nx == 0 { 0 } else { table.len() / nx.max(1) },
                actual_cols: nx,
            });
        }
        let xs_f64: Vec<f64> = xs.iter().copied().map(Into::into).collect();
        let ys_f64: Vec<f64> = ys.iter().copied().map(Into::into).collect();
        let dir_x = algo::validate_axis("x", &xs_f64)?;
        let dir_y = algo::validate_axis("y", &ys_f64)?;
        Ok(Self {
            xs,
            ys,
            table,
            nx,
            ny,
            dir_x,
            dir_y,
            y_reflect_offset: None,
            regions: Vec::new(),
            provenance: Provenance::default(),
            _markers: core::marker::PhantomData,
        })
    }

    /// Build a `Grid2D` from a uniformly-spaced **descending** y-axis and a
    /// table whose rows match that descending ordering (the common shape of
    /// publication-format tables — e.g. the Leinert (1998) zodiacal table is
    /// stored row 0 ↔ `λ−λ_sun = 180°` down to row 36 ↔ `0°`).
    ///
    /// Callers query in the same descending coordinate system in which they
    /// supplied `ys_desc`; internally the constructor:
    ///
    ///   1. Validates `ys_desc` is strictly descending **and uniform** (all
    ///      consecutive steps equal in IEEE bits).
    ///   2. Reverses the rows and the y-axis so storage is ascending.
    ///   3. Records `y_reflect_offset = ys_desc[0] + ys_desc[ny-1]`, so that
    ///      [`interp_at`](Self::interp_at) maps `y_user → y_reflect_offset −
    ///      y_user` before running the ascending-axis kernel. For uniform
    ///      grids this is bit-for-bit equivalent to having indexed the
    ///      descending axis with the legacy
    ///      `(y_max − y_user − step·i0) / step` arithmetic.
    ///
    /// The x axis is treated like any other axis (may be ascending or
    /// descending; not subject to uniformity checks).
    pub fn from_raw_row_major_y_descending(
        xs: Vec<S>,
        ys_desc: Vec<S>,
        table_desc: Vec<S>,
    ) -> Result<Self, TableError> {
        let nx = xs.len();
        let ny = ys_desc.len();
        if table_desc.len() != nx * ny {
            return Err(TableError::ShapeMismatch {
                expected_x: nx,
                expected_y: ny,
                actual_rows: if nx == 0 { 0 } else { table_desc.len() / nx.max(1) },
                actual_cols: nx,
            });
        }
        if ny < 2 {
            return Err(TableError::TooFewSamples { axis: "y", len: ny });
        }
        // Validate strictly descending + uniform (bit-equal step in IEEE 754).
        let ys_desc_f64: Vec<f64> = ys_desc.iter().copied().map(Into::into).collect();
        let step = ys_desc_f64[1] - ys_desc_f64[0];
        if !(step < 0.0) {
            return Err(TableError::NotMonotonic { axis: "y", at_index: 1 });
        }
        for i in 2..ny {
            let s = ys_desc_f64[i] - ys_desc_f64[i - 1];
            if s.to_bits() != step.to_bits() {
                return Err(TableError::NotMonotonic { axis: "y", at_index: i });
            }
        }
        // Flip y-axis to ascending storage. The reflection
        // `y_internal = (y_max + y_min) − y_user` together with reversing
        // the axis (but **keeping the rows in their original descending
        // order**) puts the storage-row index returned by ascending `locate`
        // back at the original `table_desc` row index — so no row flip is
        // required.
        let mut ys_asc: Vec<S> = ys_desc.clone();
        ys_asc.reverse();
        let table_asc: Vec<S> = table_desc;
        let xs_f64: Vec<f64> = xs.iter().copied().map(Into::into).collect();
        let dir_x = algo::validate_axis("x", &xs_f64)?;
        // Compute reflection offset = ys_desc[0] + ys_desc[ny-1] in S.
        let off = ys_desc[0] + ys_desc[ny - 1];
        Ok(Self {
            xs,
            ys: ys_asc,
            table: table_asc,
            nx,
            ny,
            dir_x,
            dir_y: AxisDirection::Ascending,
            y_reflect_offset: Some(off),
            regions: Vec::new(),
            provenance: Provenance::default(),
            _markers: core::marker::PhantomData,
        })
    }

    /// Append a [`ConstantRegion`] checked **before** interpolation. Multiple
    /// regions are checked in insertion order; the first match wins.
    pub fn with_constant_region(mut self, region: ConstantRegion<S>) -> Self {
        self.regions.push(region);
        self
    }

    /// Attach provenance metadata.
    pub fn with_provenance(mut self, provenance: Provenance) -> Self {
        self.provenance = provenance;
        self
    }

    /// Provenance metadata, if any.
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Number of x samples (columns).
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Number of y samples (rows).
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Bilinearly interpolate at the typed query point.
    ///
    /// Constant-fill regions (added via [`with_constant_region`](Self::with_constant_region))
    /// are evaluated first against the **user-supplied** `(x, y)`; if any
    /// matches, its value is returned without consulting the table. If a
    /// y-reflection is in effect (see
    /// [`from_raw_row_major_y_descending`](Self::from_raw_row_major_y_descending)),
    /// the y value passed to the bilinear kernel is `y_reflect_offset − y`.
    pub fn interp_at(
        &self,
        x: Quantity<X, S>,
        y: Quantity<Y, S>,
        oor_x: OutOfRange,
        oor_y: OutOfRange,
    ) -> Result<Quantity<V, S>, TableError> {
        let xv = x.value();
        let yv = y.value();
        for r in &self.regions {
            if r.contains(xv, yv) {
                return Ok(Quantity::<V, S>::new(r.value));
            }
        }
        let y_internal = match self.y_reflect_offset {
            Some(off) => off - yv,
            None => yv,
        };
        let xs_f64: Vec<f64> = self.xs.iter().copied().map(Into::into).collect();
        let ys_f64: Vec<f64> = self.ys.iter().copied().map(Into::into).collect();
        let table_f64: Vec<f64> = self.table.iter().copied().map(Into::into).collect();
        let rows: Vec<&[f64]> = (0..self.ny)
            .map(|iy| &table_f64[iy * self.nx..(iy + 1) * self.nx])
            .collect();
        let v = algo::bilinear(
            &xs_f64,
            &ys_f64,
            &rows,
            xv.into(),
            y_internal.into(),
            oor_x,
            oor_y,
            self.dir_x,
            self.dir_y,
        )?;
        Ok(Quantity::<V, S>::new(S::from(v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::length::{Meter, Nanometer};

    // We use Meter as a stand-in for the y-axis unit since this is just a
    // type-level test of the 2D wrapper.
    #[test]
    fn typed_corners_recover_table_values() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0];
        // Row 0: [10, 20, 30]; Row 1: [100, 200, 300]
        let table = vec![10.0, 20.0, 30.0, 100.0, 200.0, 300.0];
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major(xs, ys, table).unwrap();

        let v = g
            .interp_at(
                Quantity::<Nanometer>::new(1.0),
                Quantity::<Meter>::new(1.0),
                OutOfRange::Error,
                OutOfRange::Error,
            )
            .unwrap();
        assert_eq!(v.value(), 200.0);
    }

    #[test]
    fn rejects_shape_mismatch() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0, 1.0];
        let table = vec![1.0, 2.0, 3.0]; // wrong length
        let r: Result<Grid2D<Nanometer, Meter, Meter>, _> =
            Grid2D::from_raw_row_major(xs, ys, table);
        assert!(matches!(r, Err(TableError::ShapeMismatch { .. })));
    }

    /// Grid2D constructor accepts a descending y-axis and stores the direction.
    #[test]
    fn accepts_descending_y_axis() {
        use super::super::AxisDirection;
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![10.0_f64, 5.0, 0.0]; // descending
        let table = vec![1.0_f64; 9];
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major(xs, ys, table).unwrap();
        assert_eq!(g.dir_y, AxisDirection::Descending);
        assert_eq!(g.dir_x, AxisDirection::Ascending);
    }

    /// Critical NSB-parity test: a Grid2D with ascending x (β) and descending
    /// y (λ-λ_sun) must produce bit-for-bit identical results to NSB's manual
    /// `bilinear_unit` indexing on the same synthetic table layout.
    #[test]
    fn grid2d_descending_y_nsb_parity_bit_for_bit() {
        // Synthetic Leinert-shaped table:
        //   xs = β axis: [0, 5, 10, 15, 20] (ascending, 5 values)
        //   ys = λ-λ_sun: [20, 15, 10, 5, 0] (descending, 5 values)
        //   table[iy][ix] = (iy+1)*100 + (ix+1)*10
        let xs: Vec<f64> = vec![0.0, 5.0, 10.0, 15.0, 20.0];
        let ys: Vec<f64> = vec![20.0, 15.0, 10.0, 5.0, 0.0];
        let mut table_flat = vec![0.0_f64; 25];
        for iy in 0..5usize {
            for ix in 0..5usize {
                table_flat[iy * 5 + ix] = (iy + 1) as f64 * 100.0 + (ix + 1) as f64 * 10.0;
            }
        }

        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major(xs.clone(), ys.clone(), table_flat.clone()).unwrap();

        // Query: x=7.5 (β), y=7.5 (λ-λ_sun).
        let xq = 7.5_f64;
        let yq = 7.5_f64;

        let got = g.interp_at(
            Quantity::<Nanometer>::new(xq),
            Quantity::<Meter>::new(yq),
            OutOfRange::Error,
            OutOfRange::Error,
        ).unwrap();

        // NSB-style reference using bilinear_unit with manual descending-y indexing:
        // ix0 = floor(xq/5) = 1, bt = (7.5 - 5) / 5 = 0.5
        // ys[1]=15 >= yq=7.5 > ys[2]=10 → iy0 = 1
        // lt = (yq - ys[iy0]) / (ys[iy0+1] - ys[iy0]) = (7.5 - 15) / (10 - 15) = 1.5
        // Wait — let me recheck: ys[1]=15, ys[2]=10; yq=7.5 is NOT in [10,15], try iy0=2:
        // ys[2]=10 >= 7.5 > ys[3]=5 → iy0 = 2
        // lt = (7.5 - 10) / (5 - 10) = -2.5 / -5 = 0.5
        let ix0 = 1usize;
        let iy0 = 2usize;
        let bt = (xq - xs[ix0]) / (xs[ix0 + 1] - xs[ix0]);
        let lt = (yq - ys[iy0]) / (ys[iy0 + 1] - ys[iy0]);
        let expected = super::super::algo::bilinear_unit(
            table_flat[iy0 * 5 + ix0],
            table_flat[iy0 * 5 + ix0 + 1],
            table_flat[(iy0 + 1) * 5 + ix0],
            table_flat[(iy0 + 1) * 5 + ix0 + 1],
            bt,
            lt,
        );

        assert_eq!(
            got.value().to_bits(),
            expected.to_bits(),
            "Grid2D(descending y) must match bilinear_unit bit-for-bit: got={}, expected={}",
            got.value(),
            expected
        );
    }

    /// Constant-fill region short-circuits interpolation when the predicate
    /// matches; queries outside the region fall through to bilinear math.
    #[test]
    fn constant_region_short_circuits() {
        let xs = vec![0.0_f64, 10.0, 20.0];
        let ys = vec![0.0_f64, 10.0, 20.0];
        let table = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ];
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major(xs, ys, table)
                .unwrap()
                .with_constant_region(ConstantRegion::lower_corner(5.0, 5.0, 999.0));

        // Inside region: x<5 && y<5 → 999.
        let v_in = g.interp_at(
            Quantity::<Nanometer>::new(2.0),
            Quantity::<Meter>::new(2.0),
            OutOfRange::Error,
            OutOfRange::Error,
        ).unwrap();
        assert_eq!(v_in.value(), 999.0);

        // Boundary excluded (half-open: x<5 false at x=5).
        let v_b = g.interp_at(
            Quantity::<Nanometer>::new(5.0),
            Quantity::<Meter>::new(2.0),
            OutOfRange::Error,
            OutOfRange::Error,
        ).unwrap();
        assert_ne!(v_b.value(), 999.0);

        // Outside region: standard bilinear.
        let v_out = g.interp_at(
            Quantity::<Nanometer>::new(15.0),
            Quantity::<Meter>::new(15.0),
            OutOfRange::Error,
            OutOfRange::Error,
        ).unwrap();
        assert_eq!(v_out.value(), 7.0); // bilinear at (15,15) on this table = 7
    }

    /// Constant regions are checked in insertion order; first match wins.
    #[test]
    fn constant_regions_first_match_wins() {
        let xs = vec![0.0_f64, 10.0];
        let ys = vec![0.0_f64, 10.0];
        let table = vec![1.0, 2.0, 3.0, 4.0];
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major(xs, ys, table)
                .unwrap()
                .with_constant_region(ConstantRegion::lower_corner(5.0, 5.0, 100.0))
                .with_constant_region(ConstantRegion::lower_corner(8.0, 8.0, 200.0));

        // (2,2): both match → 100 (first).
        let v = g.interp_at(
            Quantity::<Nanometer>::new(2.0),
            Quantity::<Meter>::new(2.0),
            OutOfRange::Error, OutOfRange::Error,
        ).unwrap();
        assert_eq!(v.value(), 100.0);

        // (6,6): only second matches → 200.
        let v = g.interp_at(
            Quantity::<Nanometer>::new(6.0),
            Quantity::<Meter>::new(6.0),
            OutOfRange::Error, OutOfRange::Error,
        ).unwrap();
        assert_eq!(v.value(), 200.0);
    }

    /// `from_raw_row_major_y_descending` rejects a non-uniform descending axis.
    #[test]
    fn y_descending_rejects_non_uniform() {
        let xs = vec![0.0_f64, 1.0];
        let ys = vec![10.0_f64, 7.0, 0.0]; // steps -3 then -7 (non-uniform).
        let table = vec![1.0_f64; 6];
        let r: Result<Grid2D<Nanometer, Meter, Meter>, _> =
            Grid2D::from_raw_row_major_y_descending(xs, ys, table);
        assert!(matches!(r, Err(TableError::NotMonotonic { .. })));
    }

    /// `from_raw_row_major_y_descending` rejects an ascending axis.
    #[test]
    fn y_descending_rejects_ascending() {
        let xs = vec![0.0_f64, 1.0];
        let ys = vec![0.0_f64, 5.0, 10.0];
        let table = vec![1.0_f64; 6];
        let r: Result<Grid2D<Nanometer, Meter, Meter>, _> =
            Grid2D::from_raw_row_major_y_descending(xs, ys, table);
        assert!(matches!(r, Err(TableError::NotMonotonic { .. })));
    }

    /// y-descending constructor: querying in the descending coordinate gives
    /// the same result as the equivalent ascending-axis grid.
    #[test]
    fn y_descending_round_trip() {
        // Table laid out row-descending: row 0 ↔ y=20, row 4 ↔ y=0.
        let xs = vec![0.0_f64, 5.0, 10.0];
        let ys_desc = vec![20.0_f64, 15.0, 10.0, 5.0, 0.0];
        let mut table = Vec::with_capacity(15);
        for iy in 0..5 {
            for ix in 0..3 {
                table.push((iy * 10 + ix) as f64);
            }
        }
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major_y_descending(xs.clone(), ys_desc.clone(), table.clone())
                .unwrap();

        // Query at exact stored y values reproduces table cells.
        for (iy, &yv) in ys_desc.iter().enumerate() {
            for ix in 0..3 {
                let v = g.interp_at(
                    Quantity::<Nanometer>::new(xs[ix]),
                    Quantity::<Meter>::new(yv),
                    OutOfRange::Error, OutOfRange::Error,
                ).unwrap();
                assert_eq!(v.value(), table[iy * 3 + ix]);
            }
        }
    }

    /// y-descending bit-parity: matches the legacy NSB `(180 − dl − 5·l0)/5`
    /// arithmetic exactly. Mirrors the construction NSB performs but lets
    /// callers skip the manual `dl_asc = 180 − dl` flip.
    #[test]
    fn y_descending_matches_legacy_bit_for_bit() {
        // Synthetic Leinert-shaped table.
        let xs: Vec<f64> = (0..=18).map(|i| i as f64 * 5.0).collect();
        let ys_desc: Vec<f64> = (0..37).map(|i| 180.0 - i as f64 * 5.0).collect();
        // table[iy][ix] = iy*100 + ix*7.5 (arbitrary, distinct values).
        let mut table = Vec::with_capacity(37 * 19);
        for iy in 0..37 {
            for ix in 0..19 {
                table.push((iy as f64) * 100.0 + (ix as f64) * 7.5);
            }
        }
        let g: Grid2D<Nanometer, Meter, Meter> =
            Grid2D::from_raw_row_major_y_descending(xs.clone(), ys_desc.clone(), table.clone())
                .unwrap();

        // Pick a non-grid (β, dl) point and compare to the legacy hand-rolled
        // form: l0 = floor((180-dl)/5), lt = (180 - dl - 5·l0)/5.
        let dl = 27.3_f64;
        let beta = 12.7_f64;

        let b0 = (beta / 5.0).floor() as usize;
        let bt = (beta - 5.0 * b0 as f64) / 5.0;
        let l0_idx = ((180.0 - dl.ceil()) / 5.0).floor() as isize;
        let l0 = l0_idx.max(0).min(35) as usize;
        let lt = (180.0 - dl - 5.0 * l0 as f64) / 5.0;
        let row_lo = l0; // index in row-descending storage
        let row_hi = l0 + 1;
        let expected = super::super::algo::bilinear_unit(
            table[row_lo * 19 + b0],
            table[row_lo * 19 + b0 + 1],
            table[row_hi * 19 + b0],
            table[row_hi * 19 + b0 + 1],
            bt,
            lt,
        );

        let got = g.interp_at(
            Quantity::<Nanometer>::new(beta),
            Quantity::<Meter>::new(dl),
            OutOfRange::ClampToEndpoints,
            OutOfRange::ClampToEndpoints,
        ).unwrap();
        assert_eq!(
            got.value().to_bits(),
            expected.to_bits(),
            "y-descending Grid2D must match legacy bit-for-bit: got={}, expected={}",
            got.value(),
            expected
        );
    }
}
