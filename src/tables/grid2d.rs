// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed 2D bilinearly-interpolated table.

use crate::ext_qtty::{Quantity, Scalar, Unit};
use crate::interp::OutOfRange;
use crate::provenance::Provenance;

use super::{algo, TableError};

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
        algo::validate_axis("x", &xs_f64)?;
        algo::validate_axis("y", &ys_f64)?;
        Ok(Self {
            xs,
            ys,
            table,
            nx,
            ny,
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
}
