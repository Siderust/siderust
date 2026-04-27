// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed 1D linearly-interpolated table.

use crate::ext_qtty::{Quantity, Scalar, Unit};
use crate::interp::OutOfRange;
use crate::provenance::Provenance;

use super::{algo, TableError};

/// Strictly-monotonic 1D table of `(X, V)` samples with linear
/// interpolation. The `S` scalar parameter defaults to `f64`.
#[derive(Debug, Clone)]
pub struct Grid1D<X: Unit, V: Unit, S: Scalar = f64> {
    xs: Vec<S>,
    vs: Vec<S>,
    provenance: Provenance,
    _markers: core::marker::PhantomData<(X, V)>,
}

impl<X: Unit, V: Unit, S: Scalar + Into<f64> + From<f64>> Grid1D<X, V, S> {
    /// Build a `Grid1D` from raw scalar slices already expressed in the
    /// declared `X` / `V` units. The axis is validated to be strictly
    /// monotonic.
    pub fn from_raw(xs: Vec<S>, vs: Vec<S>) -> Result<Self, TableError> {
        if xs.len() != vs.len() {
            return Err(TableError::ShapeMismatch {
                expected_x: xs.len(),
                expected_y: 1,
                actual_rows: 1,
                actual_cols: vs.len(),
            });
        }
        let xs_f64: Vec<f64> = xs.iter().copied().map(Into::into).collect();
        algo::validate_axis("x", &xs_f64)?;
        Ok(Self {
            xs,
            vs,
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

    /// Number of samples in the table.
    pub fn len(&self) -> usize {
        self.xs.len()
    }

    /// `true` if the table has no samples.
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    /// Linearly interpolate at the typed query point.
    pub fn interp_at(
        &self,
        x: Quantity<X, S>,
        oor: OutOfRange,
    ) -> Result<Quantity<V, S>, TableError> {
        let xs_f64: Vec<f64> = self.xs.iter().copied().map(Into::into).collect();
        let vs_f64: Vec<f64> = self.vs.iter().copied().map(Into::into).collect();
        let v = algo::linear_1d(&xs_f64, &vs_f64, x.value().into(), oor)?;
        Ok(Quantity::<V, S>::new(S::from(v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::length::{Meter, Nanometer};

    #[test]
    fn typed_round_trip() {
        let xs = vec![400.0, 500.0, 600.0];
        let vs = vec![1.0, 2.0, 3.0];
        let g: Grid1D<Nanometer, Meter> = Grid1D::from_raw(xs, vs).unwrap();
        let q = Quantity::<Nanometer>::new(450.0);
        let v = g.interp_at(q, OutOfRange::Error).unwrap();
        assert_eq!(v.value(), 1.5);
    }

    #[test]
    fn rejects_non_monotonic() {
        let xs = vec![1.0_f64, 2.0, 2.0];
        let vs = vec![1.0_f64, 2.0, 3.0];
        let g: Result<Grid1D<Nanometer, Meter>, _> = Grid1D::from_raw(xs, vs);
        assert!(matches!(g, Err(TableError::NotMonotonic { .. })));
    }
}
