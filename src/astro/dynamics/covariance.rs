// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-tagged 6×6 state covariance represented as three typed 3×3 blocks.
//!
//! The covariance is decomposed into three canonical blocks in `[r, v]` order:
//!
//! ```text
//! P = [ Prr  Prv ]
//!     [ Pvr  Pvv ]
//! ```
//!
//! - `Prr`: position–position covariance (symmetric, units: km²)
//! - `Prv`: position–velocity cross-covariance (units: km·km/s)
//! - `Pvv`: velocity–velocity covariance (symmetric, units: (km/s)²)
//! - `Pvr` is derived as `Prv^T` and is never stored separately.
//!
//! Wrapped in [`StateCovariance<F>`] the blocks carry a frame tag so an
//! inertial covariance ([`StateCovariance<GCRS>`]) cannot be accidentally
//! mixed with a local-frame covariance ([`StateCovariance<RTN>`], …).
//!
//! ## Frame transport
//!
//! The Cartesian ↔ local transform uses the block-diagonal rotation
//! `T = blockdiag(R, R)`, i.e. position and velocity are rotated by the same
//! basis matrix.  In block form this reduces to three independent 3×3
//! similarity transforms (one per block), which avoids assembling the full 6×6
//! matrix.  This is the standard *instantaneous* convention and ignores the
//! time-derivative `Ṙ` of the rotation.
//!
//! ## Row-major export
//!
//! Call [`StateCovariance::to_row_major`] to obtain a `[[f64; 6]; 6]` array
//! (with `Pvr = Prv^T` filled in) for serialisation or numerical routines that
//! expect a dense 6×6 matrix.
//!
//! ## Example
//!
//! ```rust
//! use siderust::astro::dynamics::{
//!     OrbitState, Position, Velocity,
//!     covariance::StateCovariance,
//!     frames::{LocalOrbitalFrame, RTN},
//! };
//! use siderust::coordinates::frames::GCRS;
//! use siderust::time::JulianDate;
//!
//! let s = OrbitState::new(
//!     JulianDate::new(2_451_545.0),
//!     Position::<GCRS>::new(7000.0, 100.0, -200.0),
//!     Velocity::<GCRS>::new(0.5, 7.5, 0.1),
//! );
//! // Diagonal covariance: σr = 1 km, σv = 0.001 km/s.
//! let p_gcrs = StateCovariance::<GCRS>::from_stddevs([1.0, 1.0, 1.0], [1e-3, 1e-3, 1e-3]);
//! let f = LocalOrbitalFrame::<RTN>::from_state(&s);
//! let p_rtn = p_gcrs.transform_into::<RTN>(&f);
//! let p_back = p_rtn.transform_into_inertial(&f);
//! // Round-trip should recover the original covariance within floating-point error.
//! let orig = p_gcrs.to_row_major();
//! let back = p_back.to_row_major();
//! for i in 0..6 {
//!     for j in 0..6 {
//!         assert!((back[i][j] - orig[i][j]).abs() < 1e-12,
//!                 "back[{i}][{j}] = {} ≠ {}", back[i][j], orig[i][j]);
//!     }
//! }
//! ```

use affn::frames::ReferenceFrame;
use affn::matrix3::{FrameMatrix3, SymmetricFrameMatrix3};

use crate::coordinates::frames::GCRS;

use super::frames::LocalOrbitalFrame;

/// Frame-tagged 6×6 Cartesian state covariance stored as three 3×3 blocks.
///
/// The phantom parameter `F` is the frame the covariance is expressed in.
/// The three blocks are:
///
/// - `rr`: position–position (symmetric, km²)
/// - `rv`: position–velocity cross block (km·km/s)
/// - `vv`: velocity–velocity (symmetric, (km/s)²)
///
/// `vr = rv^T` is always derived on demand.
#[derive(Debug, Clone, Copy)]
pub struct StateCovariance<F> {
    rr: SymmetricFrameMatrix3<F>,
    rv: FrameMatrix3<F>,
    vv: SymmetricFrameMatrix3<F>,
}

impl<F> StateCovariance<F> {
    /// Construct from three explicit blocks.
    ///
    /// `rr` and `vv` must be symmetric (enforced by their types). `rv` is the
    /// position–velocity cross block; `vr` is derived as `rv^T`.
    #[inline]
    pub fn from_blocks(
        rr: SymmetricFrameMatrix3<F>,
        rv: FrameMatrix3<F>,
        vv: SymmetricFrameMatrix3<F>,
    ) -> Self {
        Self { rr, rv, vv }
    }

    /// Diagonal covariance built from per-axis standard deviations.
    ///
    /// `sigma_r[i]` is the 1-σ uncertainty in position component `i` (km).
    /// `sigma_v[i]` is the 1-σ uncertainty in velocity component `i` (km/s).
    /// All cross-covariances are set to zero.
    pub fn from_stddevs(sigma_r: [f64; 3], sigma_v: [f64; 3]) -> Self {
        let rr = SymmetricFrameMatrix3::from_diagonal([
            sigma_r[0] * sigma_r[0],
            sigma_r[1] * sigma_r[1],
            sigma_r[2] * sigma_r[2],
        ]);
        let vv = SymmetricFrameMatrix3::from_diagonal([
            sigma_v[0] * sigma_v[0],
            sigma_v[1] * sigma_v[1],
            sigma_v[2] * sigma_v[2],
        ]);
        Self {
            rr,
            rv: FrameMatrix3::zero(),
            vv,
        }
    }

    /// Position–position covariance block (symmetric, km²).
    #[inline]
    pub fn rr(&self) -> &SymmetricFrameMatrix3<F> {
        &self.rr
    }

    /// Position–velocity cross-covariance block (km·km/s).
    #[inline]
    pub fn rv(&self) -> &FrameMatrix3<F> {
        &self.rv
    }

    /// Velocity–position cross-covariance block (km·km/s) — derived as `rv^T`.
    #[inline]
    pub fn vr(&self) -> FrameMatrix3<F> {
        self.rv.transpose()
    }

    /// Velocity–velocity covariance block (symmetric, (km/s)²).
    #[inline]
    pub fn vv(&self) -> &SymmetricFrameMatrix3<F> {
        &self.vv
    }

    /// Assemble the full 6×6 row-major matrix `[[f64; 6]; 6]` for serialisation
    /// or numerical routines.
    ///
    /// Layout: rows/columns `0..2` are position, `3..5` are velocity.
    /// The `vr` block is filled as `rv^T`.
    pub fn to_row_major(&self) -> [[f64; 6]; 6] {
        let rr = self.rr.as_array();
        let rv = self.rv.as_array();
        let vv = self.vv.as_array();
        let mut out = [[0.0_f64; 6]; 6];
        for (i, rr_row) in rr.iter().enumerate() {
            for (j, rr_elt) in rr_row.iter().enumerate() {
                out[i][j] = *rr_elt;
                out[i][j + 3] = rv[i][j];
                out[i + 3][j] = rv[j][i]; // vr = rv^T
                out[i + 3][j + 3] = vv[i][j];
            }
        }
        out
    }

    /// Re-tag the covariance with a different frame.
    ///
    /// Unchecked: use only when the data is already expressed in `G`.
    #[inline]
    pub fn relabel<G>(self) -> StateCovariance<G> {
        StateCovariance {
            rr: self.rr.relabel::<G>(),
            rv: self.rv.relabel::<G>(),
            vv: self.vv.relabel::<G>(),
        }
    }
}

// =============================================================================
// Inertial ↔ local-frame transport
// =============================================================================

impl StateCovariance<GCRS> {
    /// Rotate this inertial covariance into the local orbital frame `M`.
    ///
    /// Uses the block-diagonal rotation `T = blockdiag(R, R)` in block form,
    /// applying `R·block·Rᵀ` independently to each of the three sub-blocks.
    pub fn transform_into<M: ReferenceFrame>(
        &self,
        frame: &LocalOrbitalFrame<M>,
    ) -> StateCovariance<M> {
        let r = frame.rotation();
        StateCovariance {
            rr: self.rr.rotated_by::<M>(&r),
            rv: self.rv.rotated_by::<M>(&r),
            vv: self.vv.rotated_by::<M>(&r),
        }
    }
}

impl<M: ReferenceFrame> StateCovariance<M> {
    /// Rotate this local-frame covariance back to the inertial parent frame
    /// ([`GCRS`]).
    pub fn transform_into_inertial(&self, frame: &LocalOrbitalFrame<M>) -> StateCovariance<GCRS> {
        let r = frame.rotation_inverse();
        StateCovariance {
            rr: self.rr.rotated_by::<GCRS>(&r),
            rv: self.rv.rotated_by::<GCRS>(&r),
            vv: self.vv.rotated_by::<GCRS>(&r),
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::frames::{LocalOrbitalFrame, RTN};
    use crate::astro::dynamics::{OrbitState, Position, Velocity};
    use crate::time::JulianDate;

    fn sample_state() -> OrbitState {
        OrbitState::new(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 100.0, -200.0),
            Velocity::<GCRS>::new(0.5, 7.5, 0.1),
        )
    }

    /// Build a physically meaningful diagonal covariance: σr = 1 km, σv = 1e-3 km/s.
    fn diagonal_covariance() -> StateCovariance<GCRS> {
        StateCovariance::from_stddevs([1.0, 1.0, 1.0], [1e-3, 1e-3, 1e-3])
    }

    #[test]
    fn from_stddevs_diagonal_layout() {
        let p = diagonal_covariance();
        // rr block: 1² = 1.0 on diagonal, 0 off-diagonal.
        assert!((p.rr().as_array()[0][0] - 1.0).abs() < 1e-15);
        assert!((p.rr().as_array()[1][1] - 1.0).abs() < 1e-15);
        assert!((p.rr().as_array()[2][2] - 1.0).abs() < 1e-15);
        assert_eq!(p.rr().as_array()[0][1], 0.0);
        // vv block: (1e-3)² = 1e-6 on diagonal.
        assert!((p.vv().as_array()[0][0] - 1e-6).abs() < 1e-21);
        // rv block: all zeros.
        for row in p.rv().as_array() {
            for v in row {
                assert_eq!(*v, 0.0);
            }
        }
    }

    #[test]
    fn vr_is_rv_transpose() {
        // Non-diagonal rv block.
        let rv_data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let p = StateCovariance::from_blocks(
            SymmetricFrameMatrix3::<GCRS>::from_diagonal([1.0, 1.0, 1.0]),
            FrameMatrix3::<GCRS>::from_array(rv_data),
            SymmetricFrameMatrix3::<GCRS>::from_diagonal([1.0, 1.0, 1.0]),
        );
        let vr = p.vr();
        for (i, row) in vr.as_array().iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert_eq!(*value, rv_data[j][i]);
            }
        }
    }

    #[test]
    fn to_row_major_block_layout() {
        let rr_data = [[2.0, 0.5, 0.0], [0.5, 3.0, 0.0], [0.0, 0.0, 1.0]];
        let rv_data = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.1]];
        let vv_data = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]];
        let p = StateCovariance::from_blocks(
            SymmetricFrameMatrix3::<GCRS>::from_upper(rr_data),
            FrameMatrix3::<GCRS>::from_array(rv_data),
            SymmetricFrameMatrix3::<GCRS>::from_upper(vv_data),
        );
        let m = p.to_row_major();
        // rr top-left
        assert_eq!(m[0][0], 2.0);
        assert_eq!(m[0][1], 0.5);
        assert_eq!(m[1][0], 0.5); // symmetric
                                  // rv top-right
        assert_eq!(m[0][3], 0.1);
        // vr bottom-left = rv^T
        assert_eq!(m[3][0], 0.1);
        // vv bottom-right
        assert_eq!(m[3][3], 0.01);
        // vr is rv^T
        for (i, row) in m[3..6].iter().enumerate() {
            for (j, value) in row[..3].iter().enumerate() {
                assert_eq!(*value, m[j][i + 3]);
            }
        }
    }

    #[test]
    fn round_trip_inertial_rtn() {
        let s = sample_state();
        let p = diagonal_covariance();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let p_rtn = p.transform_into::<RTN>(&f);
        let p_back = p_rtn.transform_into_inertial(&f);
        let orig = p.to_row_major();
        let back = p_back.to_row_major();
        for (i, row) in back.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert!(
                    (*value - orig[i][j]).abs() < 1e-12,
                    "back[{i}][{j}] = {} ≠ {}",
                    value,
                    orig[i][j]
                );
            }
        }
    }

    #[test]
    fn diagonal_covariance_stays_diagonal_under_rtn_rotation() {
        // The off-diagonal blocks of a diagonal covariance should stay at 0 under
        // block-diagonal rotation (rv block stays zero since a zero matrix rotated
        // by any matrix is still zero).
        let s = sample_state();
        let p = diagonal_covariance();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let p_rtn = p.transform_into::<RTN>(&f);
        // rv in RTN should still be zero.
        for row in p_rtn.rv().as_array() {
            for v in row {
                assert!(v.abs() < 1e-15, "rv block should remain zero: {v}");
            }
        }
    }
}
