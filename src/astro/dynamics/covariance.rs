// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-tagged 6×6 state covariance and transport between inertial and
//! local orbital frames.
//!
//! The covariance is stored as a row-major 6×6 array indexed as `[r, v]`
//! with sub-blocks
//!
//! ```text
//! [ Prr  Prv ]
//! [ Pvr  Pvv ]
//! ```
//!
//! Wrapped in [`StateCovariance<F>`] the matrix carries a frame tag in its
//! type, so an inertial covariance ([`StateCovariance<GCRS>`]) cannot be
//! accidentally mixed with a local-frame covariance
//! ([`StateCovariance<RTN>`], …).
//!
//! The Cartesian↔local transform uses the block-diagonal rotation
//! `T = blockdiag(R, R)`, i.e. position and velocity are rotated by the
//! same basis matrix.  This is the standard *instantaneous* convention and
//! ignores the time-derivative `Ṙ` of the rotation; it is valid for short
//! propagation windows.  A rate-aware variant is left for follow-up work.
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
//! let p_gcrs = StateCovariance::<GCRS>::identity();
//! let f = LocalOrbitalFrame::<RTN>::from_state(&s);
//! let p_rtn = p_gcrs.transform_into::<RTN>(&f);
//! let p_back = p_rtn.transform_into_inertial(&f);
//! for i in 0..6 {
//!     for j in 0..6 {
//!         assert!((p_back.matrix()[i][j] - p_gcrs.matrix()[i][j]).abs() < 1e-12);
//!     }
//! }
//! ```

use core::marker::PhantomData;

use affn::frames::ReferenceFrame;
use affn::Rotation3;

use crate::coordinates::frames::GCRS;

use super::frames::LocalOrbitalFrame;

/// 6×6 row-major covariance matrix (raw storage).
pub type Covariance6 = [[f64; 6]; 6];

/// Frame-tagged 6×6 Cartesian state covariance.
///
/// The phantom parameter `F` is the frame the covariance is expressed in.
#[derive(Debug, Clone, Copy)]
pub struct StateCovariance<F> {
    matrix: Covariance6,
    _frame: PhantomData<F>,
}

impl<F> StateCovariance<F> {
    /// Wrap a raw 6×6 matrix as a covariance in frame `F`.
    #[inline]
    pub fn from_matrix(matrix: Covariance6) -> Self {
        Self { matrix, _frame: PhantomData }
    }

    /// 6×6 identity covariance in frame `F`.
    pub fn identity() -> Self {
        let mut m = [[0.0; 6]; 6];
        for i in 0..6 {
            m[i][i] = 1.0;
        }
        Self::from_matrix(m)
    }

    /// Borrow the underlying row-major matrix.
    #[inline]
    pub fn matrix(&self) -> &Covariance6 {
        &self.matrix
    }

    /// Consume into the raw matrix.
    #[inline]
    pub fn into_matrix(self) -> Covariance6 {
        self.matrix
    }

    /// Re-tag the covariance with a different frame.
    ///
    /// This is an unchecked frame relabel; use sparingly and only when you
    /// know the matrix is already expressed in `G`.
    #[inline]
    pub fn relabel<G>(self) -> StateCovariance<G> {
        StateCovariance { matrix: self.matrix, _frame: PhantomData }
    }
}

// =============================================================================
// Inertial ↔ local-frame transport
// =============================================================================

impl StateCovariance<GCRS> {
    /// Rotate this inertial covariance into the local orbital frame `M`.
    pub fn transform_into<M: ReferenceFrame>(
        &self,
        frame: &LocalOrbitalFrame<M>,
    ) -> StateCovariance<M> {
        let t = block_diag(&frame.rotation());
        StateCovariance::from_matrix(similarity(&t, &self.matrix))
    }
}

impl<M: ReferenceFrame> StateCovariance<M> {
    /// Rotate this local-frame covariance back to the inertial parent
    /// frame ([`GCRS`]).
    pub fn transform_into_inertial(
        &self,
        frame: &LocalOrbitalFrame<M>,
    ) -> StateCovariance<GCRS> {
        let t = block_diag(&frame.rotation_inverse());
        StateCovariance::from_matrix(similarity(&t, &self.matrix))
    }
}

// =============================================================================
// 6×6 numeric kernels
// =============================================================================

/// Build a 6×6 block-diagonal rotation `T = blockdiag(R, R)` from a 3×3
/// rotation `R`.
///
/// The result satisfies `[r;v]_local = T · [r;v]_inertial` when `R` is the
/// inertial→local rotation.
pub fn block_diag(r: &Rotation3) -> Covariance6 {
    let m = r.as_matrix();
    let mut t = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[i][j];
            t[i + 3][j + 3] = m[i][j];
        }
    }
    t
}

/// `out = a · b`.
fn matmul6(a: &Covariance6, b: &Covariance6) -> Covariance6 {
    let mut out = [[0.0; 6]; 6];
    for i in 0..6 {
        for k in 0..6 {
            let aik = a[i][k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..6 {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

fn transpose6(a: &Covariance6) -> Covariance6 {
    let mut out = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = a[j][i];
        }
    }
    out
}

/// Rotate `P` by `T`: returns `T · P · Tᵀ`.
pub fn similarity(t: &Covariance6, p: &Covariance6) -> Covariance6 {
    let tp = matmul6(t, p);
    let tt = transpose6(t);
    matmul6(&tp, &tt)
}

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

    #[test]
    fn identity_round_trip() {
        let s = sample_state();
        let p = StateCovariance::<GCRS>::identity();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let p_rtn = p.transform_into::<RTN>(&f);
        let p_back = p_rtn.transform_into_inertial(&f);
        for i in 0..6 {
            for j in 0..6 {
                assert!((p_back.matrix()[i][j] - p.matrix()[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn identity_stays_identity_under_rotation() {
        let s = sample_state();
        let p = StateCovariance::<GCRS>::identity();
        let f = LocalOrbitalFrame::<RTN>::from_state(&s);
        let p_rtn = p.transform_into::<RTN>(&f);
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((p_rtn.matrix()[i][j] - expected).abs() < 1e-12);
            }
        }
    }
}
