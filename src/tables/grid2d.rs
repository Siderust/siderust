// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed 2D bilinearly-interpolated table.

use crate::ext_qtty::{Quantity, Scalar, Unit};
use crate::interp::OutOfRange;
use crate::provenance::Provenance;

use super::{algo, AxisDirection, TableError};

/// Strictly-monotonic 2D table of values laid out row-major as
/// `table[iy][ix]`. The two axes carry their own typed [`Unit`]
/// markers (`X` for the column axis, `Y` for the row axis); the cell
/// values carry a third [`Unit`] marker `V`.
///
/// `S` defaults to `f64`. The interpolant matches the ordering used by
/// [`algo::bilinear`] (interpolate along `x` first within each row, then
/// along `y`).
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
            provenance: Provenance::default(),
            _markers: core::marker::PhantomData,
        })
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
    pub fn interp_at(
        &self,
        x: Quantity<X, S>,
        y: Quantity<Y, S>,
        oor_x: OutOfRange,
        oor_y: OutOfRange,
    ) -> Result<Quantity<V, S>, TableError> {
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
            x.value().into(),
            y.value().into(),
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
}
