// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed 3D trilinearly-interpolated table.

use crate::ext_qtty::{Quantity, Scalar, Unit};
use crate::interp::OutOfRange;
use crate::provenance::Provenance;

use super::{algo, AxisDirection, TableError};

/// Strictly-monotonic 3D table of values with trilinear interpolation.
///
/// The three axes carry their own typed [`Unit`] markers (`X` for the
/// innermost axis, `Y` for the middle axis, `Z` for the outermost axis);
/// cell values carry a fourth [`Unit`] marker `V`.
///
/// **Storage convention** (fixed — downstream data loaders must respect this
/// when flattening their 3-D arrays into the `table` slice):
///
/// ```text
/// index = (iz · NY + iy) · NX + ix
/// ```
///
/// Innermost stride is `x`, then `y`, then `z`. Equivalently, slice
/// `table[iz * NY * NX .. (iz+1) * NY * NX]` is the z-plane at `zs[iz]`,
/// laid out row-major as `table[iy][ix]` within that plane.
///
/// `S` defaults to `f64`. The interpolant matches the ordering used by
/// [`algo::trilinear`] (interpolate along `x` first, then `y`, then `z`).
#[derive(Debug, Clone)]
pub struct Grid3D<X: Unit, Y: Unit, Z: Unit, V: Unit, S: Scalar = f64> {
    xs: Vec<S>,
    ys: Vec<S>,
    zs: Vec<S>,
    /// Flat storage; see struct doc for index convention.
    table: Vec<S>,
    nx: usize,
    ny: usize,
    nz: usize,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
    dir_z: AxisDirection,
    provenance: Provenance,
    _markers: core::marker::PhantomData<(X, Y, Z, V)>,
}

impl<X: Unit, Y: Unit, Z: Unit, V: Unit, S: Scalar + Into<f64> + From<f64>>
    Grid3D<X, Y, Z, V, S>
{
    /// Build a `Grid3D` from raw scalar slices already expressed in the
    /// declared `X` / `Y` / `Z` / `V` units. All three axes are validated
    /// to be strictly monotonic (either direction); the table is required to
    /// be `zs.len() · ys.len() · xs.len()` long, laid out with the storage
    /// convention `index = (iz·NY + iy)·NX + ix`.
    pub fn from_raw_xyz_major(
        xs: Vec<S>,
        ys: Vec<S>,
        zs: Vec<S>,
        table: Vec<S>,
    ) -> Result<Self, TableError> {
        let nx = xs.len();
        let ny = ys.len();
        let nz = zs.len();
        let expected = nx * ny * nz;
        if table.len() != expected {
            return Err(TableError::ShapeMismatch {
                expected_x: nx,
                expected_y: ny * nz,
                actual_rows: if nx == 0 { 0 } else { table.len() / nx.max(1) },
                actual_cols: nx,
            });
        }
        let xs_f64: Vec<f64> = xs.iter().copied().map(Into::into).collect();
        let ys_f64: Vec<f64> = ys.iter().copied().map(Into::into).collect();
        let zs_f64: Vec<f64> = zs.iter().copied().map(Into::into).collect();
        let dir_x = algo::validate_axis("x", &xs_f64)?;
        let dir_y = algo::validate_axis("y", &ys_f64)?;
        let dir_z = algo::validate_axis("z", &zs_f64)?;
        Ok(Self {
            xs,
            ys,
            zs,
            table,
            nx,
            ny,
            nz,
            dir_x,
            dir_y,
            dir_z,
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

    /// Number of x samples (innermost axis).
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Number of y samples (middle axis).
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Number of z samples (outermost axis).
    pub fn nz(&self) -> usize {
        self.nz
    }

    /// Trilinearly interpolate at the typed query point.
    pub fn interp_at(
        &self,
        x: Quantity<X, S>,
        y: Quantity<Y, S>,
        z: Quantity<Z, S>,
        oor_x: OutOfRange,
        oor_y: OutOfRange,
        oor_z: OutOfRange,
    ) -> Result<Quantity<V, S>, TableError> {
        let xs_f64: Vec<f64> = self.xs.iter().copied().map(Into::into).collect();
        let ys_f64: Vec<f64> = self.ys.iter().copied().map(Into::into).collect();
        let zs_f64: Vec<f64> = self.zs.iter().copied().map(Into::into).collect();
        let table_f64: Vec<f64> = self.table.iter().copied().map(Into::into).collect();
        let v = algo::trilinear(
            &xs_f64,
            &ys_f64,
            &zs_f64,
            &table_f64,
            self.nx,
            self.ny,
            x.value().into(),
            y.value().into(),
            z.value().into(),
            oor_x,
            oor_y,
            oor_z,
            self.dir_x,
            self.dir_y,
            self.dir_z,
        )?;
        Ok(Quantity::<V, S>::new(S::from(v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ext_qtty::length::{Meter, Nanometer};
    use crate::ext_qtty::time::Second;

    /// Build the canonical 3×3×3 grid where `f(x,y,z) = 100x + 10y + z`.
    fn make_333_grid() -> Grid3D<Nanometer, Meter, Second, Meter> {
        let xs: Vec<f64> = vec![0.0, 1.0, 2.0];
        let ys: Vec<f64> = vec![0.0, 1.0, 2.0];
        let zs: Vec<f64> = vec![0.0, 1.0, 2.0];
        let mut table = vec![0.0_f64; 27];
        for iz in 0..3usize {
            for iy in 0..3usize {
                for ix in 0..3usize {
                    table[(iz * 3 + iy) * 3 + ix] =
                        100.0 * ix as f64 + 10.0 * iy as f64 + iz as f64;
                }
            }
        }
        Grid3D::from_raw_xyz_major(xs, ys, zs, table).unwrap()
    }

    #[test]
    fn exact_recovery_at_integer_nodes() {
        let g = make_333_grid();
        for (xv, yv, zv) in [(0.0_f64, 0.0, 0.0), (2.0, 1.0, 0.0), (1.0, 2.0, 2.0)] {
            let expected = 100.0 * xv + 10.0 * yv + zv;
            let got = g
                .interp_at(
                    Quantity::<Nanometer>::new(xv),
                    Quantity::<Meter>::new(yv),
                    Quantity::<Second>::new(zv),
                    OutOfRange::Error,
                    OutOfRange::Error,
                    OutOfRange::Error,
                )
                .unwrap();
            assert!(
                (got.value() - expected).abs() < 1e-12,
                "at ({xv},{yv},{zv}): got {}, expected {expected}",
                got.value()
            );
        }
    }

    #[test]
    fn midpoint_gives_555() {
        // f(0.5, 0.5, 0.5) = 100*0.5 + 10*0.5 + 0.5 = 55.5
        let g = make_333_grid();
        let v = g
            .interp_at(
                Quantity::<Nanometer>::new(0.5),
                Quantity::<Meter>::new(0.5),
                Quantity::<Second>::new(0.5),
                OutOfRange::Error,
                OutOfRange::Error,
                OutOfRange::Error,
            )
            .unwrap();
        assert!((v.value() - 55.5).abs() < 1e-12, "got {}", v.value());
    }

    #[test]
    fn clamp_policy_per_axis() {
        let g = make_333_grid();
        // x OOR (above), clamp → same as x=2.0
        let v_clamp = g
            .interp_at(
                Quantity::<Nanometer>::new(5.0),
                Quantity::<Meter>::new(1.0),
                Quantity::<Second>::new(1.0),
                OutOfRange::ClampToEndpoints,
                OutOfRange::Error,
                OutOfRange::Error,
            )
            .unwrap();
        let v_exact = g
            .interp_at(
                Quantity::<Nanometer>::new(2.0),
                Quantity::<Meter>::new(1.0),
                Quantity::<Second>::new(1.0),
                OutOfRange::Error,
                OutOfRange::Error,
                OutOfRange::Error,
            )
            .unwrap();
        assert_eq!(v_clamp.value().to_bits(), v_exact.value().to_bits());
    }

    #[test]
    fn zero_policy_per_axis() {
        let g = make_333_grid();
        // z OOR (below), Zero policy → 0
        let v = g
            .interp_at(
                Quantity::<Nanometer>::new(1.0),
                Quantity::<Meter>::new(1.0),
                Quantity::<Second>::new(-1.0),
                OutOfRange::Error,
                OutOfRange::Error,
                OutOfRange::Zero,
            )
            .unwrap();
        assert_eq!(v.value(), 0.0);
    }

    #[test]
    fn error_policy_returns_out_of_range() {
        let g = make_333_grid();
        // y OOR
        let r = g.interp_at(
            Quantity::<Nanometer>::new(1.0),
            Quantity::<Meter>::new(5.0),
            Quantity::<Second>::new(1.0),
            OutOfRange::Error,
            OutOfRange::Error,
            OutOfRange::Error,
        );
        assert!(
            matches!(r, Err(TableError::OutOfRange { axis: "y", .. })),
            "expected OutOfRange on y, got {r:?}"
        );
        // z OOR
        let r = g.interp_at(
            Quantity::<Nanometer>::new(1.0),
            Quantity::<Meter>::new(1.0),
            Quantity::<Second>::new(10.0),
            OutOfRange::Error,
            OutOfRange::Error,
            OutOfRange::Error,
        );
        assert!(
            matches!(r, Err(TableError::OutOfRange { axis: "z", .. })),
            "expected OutOfRange on z, got {r:?}"
        );
    }

    #[test]
    fn mixed_oor_policies() {
        // x: Zero, y: ClampToEndpoints, z: Error — all in one call
        let g = make_333_grid();
        let v = g
            .interp_at(
                Quantity::<Nanometer>::new(-1.0), // x OOR → Zero wins
                Quantity::<Meter>::new(5.0),      // y OOR but clamp (irrelevant, zero wins)
                Quantity::<Second>::new(1.0),
                OutOfRange::Zero,
                OutOfRange::ClampToEndpoints,
                OutOfRange::Error,
            )
            .unwrap();
        assert_eq!(v.value(), 0.0);
    }

    #[test]
    fn rejects_shape_mismatch() {
        let xs = vec![0.0_f64, 1.0];
        let ys = vec![0.0_f64, 1.0];
        let zs = vec![0.0_f64, 1.0];
        let table = vec![1.0_f64; 5]; // wrong: should be 2*2*2 = 8
        let r: Result<Grid3D<Nanometer, Meter, Second, Meter>, _> =
            Grid3D::from_raw_xyz_major(xs, ys, zs, table);
        assert!(matches!(r, Err(TableError::ShapeMismatch { .. })));
    }

    #[test]
    fn typed_interp_returns_correct_unit() {
        // Verify that three distinct Unit markers compile and the return type
        // carries the V marker.
        let g = make_333_grid();
        let result: Quantity<Meter> = g
            .interp_at(
                Quantity::<Nanometer>::new(1.0),
                Quantity::<Meter>::new(1.0),
                Quantity::<Second>::new(1.0),
                OutOfRange::Error,
                OutOfRange::Error,
                OutOfRange::Error,
            )
            .unwrap();
        // f(1,1,1) = 100 + 10 + 1 = 111
        assert!((result.value() - 111.0).abs() < 1e-12);
    }
}
