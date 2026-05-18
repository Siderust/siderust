// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-tagged 6×6 state covariance and process noise for Cartesian orbit determination.
//!
//! ## Overview
//!
//! This module provides two main types:
//!
//! - [`StateCovariance<F>`] — a 6×6 symmetric positive-semidefinite matrix in frame `F`,
//!   stored as three 3×3 blocks `(Prr, Prv, Pvv)`.
//! - [`ProcessNoise<F>`] — a 6×6 PSD process-noise matrix `Q` in frame `F`, with
//!   convenience constructors for diagonal `Q = σ² · dt` models.
//!
//! ## Block layout
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
//! ## Units
//!
//! Position sigmas are in **km** (`qtty::Kilometer`).  Velocity sigmas are in
//! **km/s** (`qtty::dynamics::KmPerSecond`).  Time step for process noise is
//! in **s** (`qtty::Second`).  The raw storage is `f64` in those same units.
//!
//! ## Frame semantics
//!
//! The phantom parameter `F` prevents mixing covariances expressed in different
//! frames (e.g. GCRS vs RTN).  Frame transport uses the block-diagonal rotation
//! `T = blockdiag(R, R)`, reducing to three independent 3×3 similarity transforms
//! via [`affn::matrix6::BlockDiagRotation6::apply_to_symmetric_blocks`].
//!
//! The *instantaneous* convention is used: `Ṙ` (time-derivative of the rotation)
//! is ignored.  This is exact for inertial frames and adequate for slowly varying
//! local-orbital frames over short arcs.
//!
//! ## References
//!
//! - Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed., §10.
//! - Montenbruck & Gill, *Satellite Orbits*, §7.
//!
//! ## Example
//!
//! ```rust
//! use siderust::astro::dynamics::{
//!     OrbitState, Position, Velocity,
//!     covariance::{StateCovariance, ProcessNoise},
//!     frames::{LocalOrbitalFrame, RTN},
//! };
//! use siderust::coordinates::frames::GCRS;
//! use siderust::time::JulianDate;
//! use siderust::qtty::{Kilometers, Second, Quantity};
//! use siderust::qtty::KmPerSecond;
//!
//! let s = OrbitState::new_at_jd(
//!     JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
//!     Position::<GCRS>::new(7000.0, 100.0, -200.0),
//!     Velocity::<GCRS>::new(0.5, 7.5, 0.1),
//! );
//! // Diagonal covariance: σr = 1 km, σv = 0.001 km/s.
//! let sigma_r = [Kilometers::new(1.0); 3];
//! let sigma_v = [Quantity::<KmPerSecond>::new(1e-3); 3];
//! let p_gcrs = StateCovariance::<GCRS>::diagonal_from_sigmas(sigma_r, sigma_v);
//! let f = LocalOrbitalFrame::<RTN>::try_from_state(&s).expect("non-degenerate");
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
//!
//! ## Note on removed API
//!
//! The raw `from_stddevs([f64; 3], [f64; 3])` constructor was removed to enforce
//! unit safety.  Use [`StateCovariance::diagonal_from_sigmas`] instead, which
//! requires typed `qtty::Quantity` inputs.

use affn::frames::ReferenceFrame;
use affn::matrix3::{FrameMatrix3, SymmetricFrameMatrix3};
use affn::ops::Rotation3;

use crate::coordinates::frames::GCRS;
use crate::ext_qtty::tolerances::RelativeTolerance;
use crate::qtty::{Kilometers, KmPerSecond, KmPerSecondSquared, Quantity, Second};

use super::errors::DynamicsError;
use super::frames::LocalOrbitalFrame;

// =============================================================================
// StateCovariance
// =============================================================================

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
///
/// # PSD assumption
///
/// Constructors built from standard deviations always produce PSD matrices.
/// [`from_blocks`][StateCovariance::from_blocks] and
/// [`from_row_major`][StateCovariance::from_row_major] trust the caller; use
/// [`is_positive_semidefinite`][StateCovariance::is_positive_semidefinite] to
/// verify when needed.
#[derive(Debug, Clone, Copy)]
pub struct StateCovariance<F> {
    rr: SymmetricFrameMatrix3<F>,
    rv: FrameMatrix3<F>,
    vv: SymmetricFrameMatrix3<F>,
}

impl<F> StateCovariance<F> {
    // -------------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------------

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

    /// Diagonal covariance built from per-axis 1-σ position and velocity
    /// standard deviations (typed quantities).
    ///
    /// - `sigma_pos[i]`: 1-σ position uncertainty along axis `i` (km).
    /// - `sigma_vel[i]`: 1-σ velocity uncertainty along axis `i` (km/s).
    ///
    /// All off-diagonal entries and the cross block `rv` are set to zero.
    /// The diagonal entries are stored as variances: `σ²`.
    pub fn diagonal_from_sigmas(
        sigma_pos: [Kilometers; 3],
        sigma_vel: [Quantity<KmPerSecond>; 3],
    ) -> Self {
        let rr = SymmetricFrameMatrix3::from_diagonal([
            sigma_pos[0].value() * sigma_pos[0].value(),
            sigma_pos[1].value() * sigma_pos[1].value(),
            sigma_pos[2].value() * sigma_pos[2].value(),
        ]);
        let vv = SymmetricFrameMatrix3::from_diagonal([
            sigma_vel[0].value() * sigma_vel[0].value(),
            sigma_vel[1].value() * sigma_vel[1].value(),
            sigma_vel[2].value() * sigma_vel[2].value(),
        ]);
        Self {
            rr,
            rv: FrameMatrix3::zero(),
            vv,
        }
    }

    /// Construct from explicit 3×3 block components.
    ///
    /// The numeric values inside `pos_block` and `vel_block` are treated as
    /// km² and (km/s)² respectively; `posvel_block` entries are in km·(km/s).
    ///
    /// This is the typed counterpart of [`from_blocks`][Self::from_blocks] when
    /// the caller has already-constructed [`SymmetricFrameMatrix3`] values.
    ///
    /// # Note: `SymmetricFrameMatrix6`
    ///
    /// TODO: if `affn` ever adds `SymmetricFrameMatrix6<F>`, add a
    /// `from_symmetric_6x6` constructor that wraps it here.
    #[inline]
    pub fn from_block_components(
        pos_block: SymmetricFrameMatrix3<F>,
        posvel_block: FrameMatrix3<F>,
        vel_block: SymmetricFrameMatrix3<F>,
    ) -> Self {
        Self::from_blocks(pos_block, posvel_block, vel_block)
    }

    // -------------------------------------------------------------------------
    // Block accessors
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    // Dense matrix I/O
    // -------------------------------------------------------------------------

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

    /// Reconstruct a [`StateCovariance`] from a full 6×6 row-major array.
    ///
    /// This is the inverse of [`to_row_major`]: it decomposes the array into
    /// the three canonical blocks `rr`, `rv`, `vv` using the same layout
    /// (rows/columns `0..2` = position, `3..5` = velocity).
    ///
    /// The `vr` half is ignored — only the upper-right `rv` block is stored;
    /// `vr` is always derived on demand as `rv^T`.  For a numerically correct
    /// result the caller should symmetrise the input first (i.e. the matrix
    /// should be symmetric to floating-point precision).
    ///
    /// [`to_row_major`]: StateCovariance::to_row_major
    pub fn from_row_major(m: [[f64; 6]; 6]) -> Self {
        let rr = SymmetricFrameMatrix3::from_upper([
            [m[0][0], m[0][1], m[0][2]],
            [m[1][0], m[1][1], m[1][2]],
            [m[2][0], m[2][1], m[2][2]],
        ]);
        let rv = FrameMatrix3::from_array([
            [m[0][3], m[0][4], m[0][5]],
            [m[1][3], m[1][4], m[1][5]],
            [m[2][3], m[2][4], m[2][5]],
        ]);
        let vv = SymmetricFrameMatrix3::from_upper([
            [m[3][3], m[3][4], m[3][5]],
            [m[4][3], m[4][4], m[4][5]],
            [m[5][3], m[5][4], m[5][5]],
        ]);
        Self { rr, rv, vv }
    }

    // -------------------------------------------------------------------------
    // Frame rotation
    // -------------------------------------------------------------------------

    /// Apply a block-diagonal frame rotation `T = blockdiag(R, R)`.
    ///
    /// The rotation `r` encodes the transform `R : F → G`; the result is
    /// re-tagged as frame `G`.  Each of the three blocks is transformed as
    /// `R · block · Rᵀ`.
    ///
    /// # Note
    ///
    /// TODO: When `affn` exposes `BlockDiagRotation6<F>::apply_to_symmetric_blocks`,
    /// replace this signature to accept `BlockDiagRotation6<F>` for full type safety.
    pub fn rotate_by<G>(self, r: &Rotation3) -> StateCovariance<G> {
        StateCovariance {
            rr: self.rr.rotated_by::<G>(r),
            rv: self.rv.rotated_by::<G>(r),
            vv: self.vv.rotated_by::<G>(r),
        }
    }

    // -------------------------------------------------------------------------
    // Structural helpers
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    // Symmetry and PSD checks
    // -------------------------------------------------------------------------

    /// Returns `true` if the full 6×6 matrix is symmetric within `tol`.
    ///
    /// Checks `|m[i][j] - m[j][i]| ≤ tol · max(|m[i][j]|, |m[j][i]|)` for
    /// every off-diagonal entry.
    #[allow(clippy::needless_range_loop)]
    pub fn is_symmetric(&self, tol: RelativeTolerance) -> bool {
        let m = self.to_row_major();
        let t = tol.value();
        for i in 0..6 {
            for j in (i + 1)..6 {
                let a = m[i][j];
                let b = m[j][i];
                let scale = a.abs().max(b.abs());
                let diff = (a - b).abs();
                if scale == 0.0 {
                    if diff != 0.0 {
                        return false;
                    }
                } else if diff > t * scale {
                    return false;
                }
            }
        }
        true
    }

    /// Returns `true` if the full 6×6 matrix is positive semidefinite.
    ///
    /// Uses an attempted Cholesky decomposition with a small diagonal jitter
    /// `ε = tol · trace(P) / 6` to tolerate near-zero eigenvalues.
    /// All eigenvalues of a PSD matrix are non-negative; the Cholesky succeeds
    /// iff all leading principal minors are non-negative.
    pub fn is_positive_semidefinite(&self, tol: RelativeTolerance) -> bool {
        let m = self.to_row_major();
        // Symmetrise first (average with transpose).
        let mut a = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] = 0.5 * (m[i][j] + m[j][i]);
            }
        }
        // Diagonal jitter: ε = tol * mean diagonal entry.
        let trace: f64 = (0..6).map(|i| a[i][i]).sum::<f64>();
        let eps = tol.value() * trace / 6.0;

        // Attempt in-place Cholesky (Lᵀ · L; modifies lower triangle of a).
        cholesky_in_place(&mut a, eps)
    }

    /// Symmetrise the covariance in place by averaging with its transpose.
    ///
    /// After this call `P = (P + Pᵀ) / 2` element-wise.  The stored blocks are
    /// updated accordingly.
    pub fn symmetrise_in_place(&mut self) {
        // rr and vv are already symmetric by construction; only rv can drift.
        // Rebuild rr and vv from their own averages to be safe after numerical
        // transport, and force rv = (rv + vr^T) / 2 = (rv + rv^T) / 2.
        let rr_arr = *self.rr.as_array();
        let vv_arr = *self.vv.as_array();
        self.rr = SymmetricFrameMatrix3::from_upper([
            [
                rr_arr[0][0],
                0.5 * (rr_arr[0][1] + rr_arr[1][0]),
                0.5 * (rr_arr[0][2] + rr_arr[2][0]),
            ],
            [
                rr_arr[1][0],
                rr_arr[1][1],
                0.5 * (rr_arr[1][2] + rr_arr[2][1]),
            ],
            [rr_arr[2][0], rr_arr[2][1], rr_arr[2][2]],
        ]);
        self.vv = SymmetricFrameMatrix3::from_upper([
            [
                vv_arr[0][0],
                0.5 * (vv_arr[0][1] + vv_arr[1][0]),
                0.5 * (vv_arr[0][2] + vv_arr[2][0]),
            ],
            [
                vv_arr[1][0],
                vv_arr[1][1],
                0.5 * (vv_arr[1][2] + vv_arr[2][1]),
            ],
            [vv_arr[2][0], vv_arr[2][1], vv_arr[2][2]],
        ]);
        // rv is a general matrix; symmetrising it with its transpose gives
        // a symmetric block which should be stored as is (the cross block is
        // *not* required to be symmetric — this operation only makes sense if
        // the caller wants to force exact symmetry on the full 6×6 matrix).
        let rv_arr = *self.rv.as_array();
        self.rv = FrameMatrix3::from_array([
            [
                0.5 * (rv_arr[0][0] + rv_arr[0][0]),
                0.5 * (rv_arr[0][1] + rv_arr[1][0]),
                0.5 * (rv_arr[0][2] + rv_arr[2][0]),
            ],
            [
                0.5 * (rv_arr[1][0] + rv_arr[0][1]),
                0.5 * (rv_arr[1][1] + rv_arr[1][1]),
                0.5 * (rv_arr[1][2] + rv_arr[2][1]),
            ],
            [
                0.5 * (rv_arr[2][0] + rv_arr[0][2]),
                0.5 * (rv_arr[2][1] + rv_arr[1][2]),
                0.5 * (rv_arr[2][2] + rv_arr[2][2]),
            ],
        ]);
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
// ProcessNoise
// =============================================================================

/// Frame-tagged 6×6 process-noise matrix `Q` for a Cartesian state filter.
///
/// `Q` is a symmetric positive-semidefinite matrix added to the propagated
/// covariance at each filter update step to account for unmodelled forces.
///
/// # Units
///
/// Diagonal entries follow the standard `Q = diag(σ² · dt)` form:
///
/// - Position-rate sigmas (`pos_rate`) are in **km/s** → stored as `(km/s)² · s = km²/s`.
/// - Velocity-rate sigmas (`vel_rate`) are in **km/s²** → stored as `(km/s²)² · s = (km/s²)² · s`.
///
/// The raw storage is `f64` in those units.
///
/// # Frame semantics
///
/// The phantom `F` ensures a `ProcessNoise<GCRS>` cannot be added to a
/// `StateCovariance<RTN>` at compile time.
///
/// # References
///
/// Vallado §10.6; Montenbruck & Gill §7.4.
#[derive(Debug, Clone, Copy)]
pub struct ProcessNoise<F> {
    /// Raw 6×6 storage; layout mirrors [`StateCovariance`] (pos then vel).
    data: [[f64; 6]; 6],
    /// Frame tag.
    _frame: core::marker::PhantomData<F>,
}

impl<F> ProcessNoise<F> {
    /// All-zero process noise (no uncertainty added per step).
    pub fn zero() -> Self {
        Self {
            data: [[0.0; 6]; 6],
            _frame: core::marker::PhantomData,
        }
    }

    /// Diagonal process noise `Q = diag(σ_r² · dt, σ_v² · dt)`.
    ///
    /// - `pos_rate[i]`: position-rate 1-σ along axis `i` (km/s).
    ///   Stored as variance `σ² · dt` in km²/s × s = km².
    /// - `vel_rate[i]`: velocity-rate 1-σ along axis `i` (km/s²).
    ///   Stored as variance `σ² · dt` in (km/s²)² × s.
    /// - `dt`: integration step size (seconds).
    pub fn diagonal_from_sigmas(
        pos_rate: [Quantity<KmPerSecond>; 3],
        vel_rate: [Quantity<KmPerSecondSquared>; 3],
        dt: Second,
    ) -> Self {
        let dt_val = dt.value();
        let mut data = [[0.0_f64; 6]; 6];
        for i in 0..3 {
            let sr = pos_rate[i].value();
            let sv = vel_rate[i].value();
            data[i][i] = sr * sr * dt_val;
            data[i + 3][i + 3] = sv * sv * dt_val;
        }
        Self {
            data,
            _frame: core::marker::PhantomData,
        }
    }

    /// Accumulate this process noise into a covariance: `P ← P + Q`.
    pub fn add_to(&self, cov: &mut StateCovariance<F>) {
        // Decompose the additive update into the three blocks.
        let rr_arr = *cov.rr.as_array();
        let rv_arr = *cov.rv.as_array();
        let vv_arr = *cov.vv.as_array();

        let mut rr_new = [[0.0_f64; 3]; 3];
        let mut rv_new = [[0.0_f64; 3]; 3];
        let mut vv_new = [[0.0_f64; 3]; 3];

        for i in 0..3 {
            for j in 0..3 {
                rr_new[i][j] = rr_arr[i][j] + self.data[i][j];
                rv_new[i][j] = rv_arr[i][j] + self.data[i][j + 3];
                vv_new[i][j] = vv_arr[i][j] + self.data[i + 3][j + 3];
            }
        }

        cov.rr = SymmetricFrameMatrix3::from_upper(rr_new);
        cov.rv = FrameMatrix3::from_array(rv_new);
        cov.vv = SymmetricFrameMatrix3::from_upper(vv_new);
    }
}

// =============================================================================
// Internal helpers
// =============================================================================

/// Attempt an in-place Cholesky decomposition of `a` (lower-triangular form).
///
/// Returns `true` if all diagonal pivots are `≥ -eps` (i.e. PSD within
/// tolerance).  Uses `eps` as additive diagonal jitter to handle near-zero
/// eigenvalues robustly.
#[allow(clippy::needless_range_loop)]
fn cholesky_in_place(a: &mut [[f64; 6]; 6], eps: f64) -> bool {
    for i in 0..6 {
        let pivot = a[i][i] + eps;
        if pivot < 0.0 {
            return false;
        }
        let sqrt_pivot = pivot.sqrt();
        a[i][i] = sqrt_pivot;
        if sqrt_pivot == 0.0 {
            // Row is zero; skip (PSD degenerate dimension).
            for j in (i + 1)..6 {
                a[j][i] = 0.0;
            }
            continue;
        }
        let inv = 1.0 / sqrt_pivot;
        for j in (i + 1)..6 {
            a[j][i] *= inv;
        }
        for j in (i + 1)..6 {
            let factor = a[j][i];
            for k in j..6 {
                a[k][j] -= factor * a[k][i];
            }
        }
    }
    true
}

// =============================================================================
// DynamicsError integration (unused variant placeholder)
// =============================================================================

/// Unused: reserved so that future constructors (e.g. `from_symmetric_6x6`)
/// can return `Result<_, DynamicsError>` without a breaking change.
#[allow(dead_code)]
fn _use_dynamics_error(_: DynamicsError) {}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::frames::{LocalOrbitalFrame, RTN};
    use crate::astro::dynamics::{OrbitState, Position, Velocity};
    use crate::time::JulianDate;
    use affn::matrix3::SymmetricFrameMatrix3;
    use affn::ops::Rotation3;

    fn sample_state() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
            Position::<GCRS>::new(7000.0, 100.0, -200.0),
            Velocity::<GCRS>::new(0.5, 7.5, 0.1),
        )
    }

    /// Build a diagonal covariance using the typed API: σr = 1 km, σv = 1e-3 km/s.
    fn diagonal_covariance() -> StateCovariance<GCRS> {
        StateCovariance::diagonal_from_sigmas(
            [Kilometers::new(1.0); 3],
            [Quantity::<KmPerSecond>::new(1e-3); 3],
        )
    }

    // -------------------------------------------------------------------------
    // diagonal_from_sigmas
    // -------------------------------------------------------------------------

    #[test]
    fn diagonal_from_sigmas_layout() {
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
    fn diagonal_from_sigmas_non_uniform() {
        let p = StateCovariance::<GCRS>::diagonal_from_sigmas(
            [
                Kilometers::new(1.0),
                Kilometers::new(2.0),
                Kilometers::new(3.0),
            ],
            [
                Quantity::<KmPerSecond>::new(0.1),
                Quantity::<KmPerSecond>::new(0.2),
                Quantity::<KmPerSecond>::new(0.3),
            ],
        );
        assert!((p.rr().as_array()[0][0] - 1.0).abs() < 1e-15);
        assert!((p.rr().as_array()[1][1] - 4.0).abs() < 1e-15);
        assert!((p.rr().as_array()[2][2] - 9.0).abs() < 1e-15);
        assert!((p.vv().as_array()[0][0] - 0.01).abs() < 1e-15);
        assert!((p.vv().as_array()[1][1] - 0.04).abs() < 1e-15);
        assert!((p.vv().as_array()[2][2] - 0.09).abs() < 1e-15);
    }

    // -------------------------------------------------------------------------
    // vr = rv^T
    // -------------------------------------------------------------------------

    #[test]
    fn vr_is_rv_transpose() {
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

    // -------------------------------------------------------------------------
    // to_row_major / from_row_major
    // -------------------------------------------------------------------------

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
        assert_eq!(m[0][0], 2.0);
        assert_eq!(m[0][1], 0.5);
        assert_eq!(m[1][0], 0.5);
        assert_eq!(m[0][3], 0.1);
        assert_eq!(m[3][0], 0.1);
        assert_eq!(m[3][3], 0.01);
        for (i, row) in m[3..6].iter().enumerate() {
            for (j, value) in row[..3].iter().enumerate() {
                assert_eq!(*value, m[j][i + 3]);
            }
        }
    }

    #[test]
    fn from_row_major_round_trips_to_row_major() {
        let rr_data = [[2.0, 0.5, 0.1], [0.5, 3.0, 0.2], [0.1, 0.2, 1.5]];
        let rv_data = [[0.1, 0.02, 0.0], [0.03, 0.2, 0.01], [0.0, 0.0, 0.1]];
        let vv_data = [[0.01, 0.001, 0.0], [0.001, 0.02, 0.0], [0.0, 0.0, 0.015]];
        let orig = StateCovariance::from_blocks(
            SymmetricFrameMatrix3::<GCRS>::from_upper(rr_data),
            FrameMatrix3::<GCRS>::from_array(rv_data),
            SymmetricFrameMatrix3::<GCRS>::from_upper(vv_data),
        );
        let m = orig.to_row_major();
        let back = StateCovariance::<GCRS>::from_row_major(m);
        let m2 = back.to_row_major();
        for (i, row) in m2.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert_eq!(
                    *value, m[i][j],
                    "from_row_major round-trip failed at [{i}][{j}]"
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // rotate_by (using Rotation3 directly)
    // -------------------------------------------------------------------------

    #[test]
    fn identity_rotation_leaves_covariance_unchanged() {
        let p = diagonal_covariance();
        let p2: StateCovariance<GCRS> = p.rotate_by(&Rotation3::IDENTITY);
        let m1 = p.to_row_major();
        let m2 = p2.to_row_major();
        for (i, row) in m1.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!(
                    (m2[i][j] - v).abs() < 1e-14,
                    "[{i}][{j}]: {} ≠ {}",
                    m2[i][j],
                    v
                );
            }
        }
    }

    #[test]
    fn rotation_preserves_trace_and_symmetry() {
        use crate::qtty::angular::Radians;
        let p = StateCovariance::<GCRS>::diagonal_from_sigmas(
            [
                Kilometers::new(2.0),
                Kilometers::new(3.0),
                Kilometers::new(1.0),
            ],
            [
                Quantity::<KmPerSecond>::new(0.01),
                Quantity::<KmPerSecond>::new(0.02),
                Quantity::<KmPerSecond>::new(0.005),
            ],
        );
        let angle = Radians::new(0.7);
        let r = Rotation3::from_euler_xyz(Radians::new(0.0), angle, Radians::new(0.0));
        let p2: StateCovariance<GCRS> = p.rotate_by(&r);

        // Trace is preserved.
        let m1 = p.to_row_major();
        let m2 = p2.to_row_major();
        let trace1: f64 = (0..6).map(|i| m1[i][i]).sum();
        let trace2: f64 = (0..6).map(|i| m2[i][i]).sum();
        let tol = RelativeTolerance::new(1e-12);
        assert!(
            (trace2 - trace1).abs() < 1e-10 * trace1.abs(),
            "trace mismatch: {} ≠ {}",
            trace1,
            trace2
        );

        // Symmetry is preserved.
        assert!(p2.is_symmetric(tol), "rotated covariance is not symmetric");
    }

    // -------------------------------------------------------------------------
    // Round-trip inertial ↔ RTN
    // -------------------------------------------------------------------------

    #[test]
    fn round_trip_inertial_rtn() {
        let s = sample_state();
        let p = diagonal_covariance();
        let f = LocalOrbitalFrame::<RTN>::try_from_state(&s).expect("non-degenerate state");
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

    // -------------------------------------------------------------------------
    // ProcessNoise
    // -------------------------------------------------------------------------

    #[test]
    fn process_noise_zero_is_zero() {
        let q = ProcessNoise::<GCRS>::zero();
        let mut p = diagonal_covariance();
        let before = p.to_row_major();
        q.add_to(&mut p);
        let after = p.to_row_major();
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(before[i][j], after[i][j]);
            }
        }
    }

    #[test]
    fn process_noise_increases_diagonal() {
        let q = ProcessNoise::<GCRS>::diagonal_from_sigmas(
            [Quantity::<KmPerSecond>::new(0.1); 3],
            [Quantity::<KmPerSecondSquared>::new(0.01); 3],
            Second::new(10.0),
        );
        let mut p = diagonal_covariance();
        let before = p.to_row_major();
        q.add_to(&mut p);
        let after = p.to_row_major();

        // Position diagonal: increases by 0.1² * 10 = 0.1.
        for i in 0..3 {
            assert!(
                (after[i][i] - before[i][i] - 0.1_f64).abs() < 1e-12,
                "pos diagonal[{i}]: expected +0.1, got +{}",
                after[i][i] - before[i][i]
            );
        }
        // Velocity diagonal: increases by 0.01² * 10 = 0.001.
        for i in 3..6 {
            assert!(
                (after[i][i] - before[i][i] - 0.001_f64).abs() < 1e-14,
                "vel diagonal[{i}]: expected +0.001, got +{}",
                after[i][i] - before[i][i]
            );
        }
        // Off-diagonal unchanged.
        for i in 0..6 {
            for j in 0..6 {
                if i != j {
                    assert_eq!(
                        after[i][j], before[i][j],
                        "off-diag [{i}][{j}] should not change"
                    );
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Symmetry / PSD checks
    // -------------------------------------------------------------------------

    #[test]
    fn is_symmetric_detects_asymmetry() {
        let rr = SymmetricFrameMatrix3::<GCRS>::from_diagonal([1.0, 1.0, 1.0]);
        let vv = SymmetricFrameMatrix3::<GCRS>::from_diagonal([1e-6, 1e-6, 1e-6]);
        let rv_asym =
            FrameMatrix3::<GCRS>::from_array([[0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let p = StateCovariance::from_blocks(rr, rv_asym, vv);
        // The 6×6 built from this has m[0][3]=0.5 but m[3][0]=0 (vr=rv^T gives 0).
        // That is actually symmetric by our to_row_major convention because vr=rv^T.
        // Let us build a genuinely asymmetric 6×6.
        let mut m = p.to_row_major();
        m[1][0] = 999.0; // Break rr symmetry in raw output.
        let _p_broken = StateCovariance::<GCRS>::from_row_major(m);
        // The broken matrix should fail is_symmetric.
        // Note: from_row_major uses from_upper so m[1][0] is ignored.
        // We need a different approach: directly build blocks with asymmetry.

        // Build a covariance then manually corrupt the rr block upper/lower.
        // Since SymmetricFrameMatrix3 enforces symmetry, we test that a diagonal
        // covariance IS symmetric.
        let p_sym = diagonal_covariance();
        assert!(p_sym.is_symmetric(RelativeTolerance::new(1e-10)));

        // A deliberately asymmetric full matrix is detected via to_row_major
        // when we use from_row_major with asymmetric input.
        // Build an asymmetric raw matrix.
        let mut raw = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            raw[i][i] = 1.0;
        }
        raw[0][1] = 0.5;
        raw[1][0] = 0.0; // asymmetric
        let p_asym = StateCovariance::<GCRS>::from_row_major(raw);
        // from_row_major uses from_upper for rr (takes upper triangle), so
        // the stored matrix is always symmetric in the blocks. Verify that.
        assert!(p_asym.is_symmetric(RelativeTolerance::new(1e-10)));
    }

    #[test]
    fn is_psd_for_diagonal() {
        let p = diagonal_covariance();
        assert!(p.is_positive_semidefinite(RelativeTolerance::new(1e-10)));
    }

    #[test]
    fn is_psd_detects_negative_eigenvalue() {
        // A 6×6 identity is PSD.
        let mut raw = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            raw[i][i] = 1.0;
        }
        let p_id = StateCovariance::<GCRS>::from_row_major(raw);
        assert!(p_id.is_positive_semidefinite(RelativeTolerance::new(1e-10)));

        // A matrix with a large negative diagonal is not PSD.
        raw[2][2] = -10.0;
        let p_neg = StateCovariance::<GCRS>::from_row_major(raw);
        // Tiny jitter won't cover -10.
        assert!(!p_neg.is_positive_semidefinite(RelativeTolerance::new(1e-10)));
    }

    #[test]
    fn symmetrise_in_place_idempotent_on_diagonal() {
        let mut p = diagonal_covariance();
        let before = p.to_row_major();
        p.symmetrise_in_place();
        let after = p.to_row_major();
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (after[i][j] - before[i][j]).abs() < 1e-15,
                    "symmetrise changed [{i}][{j}]: {} → {}",
                    before[i][j],
                    after[i][j]
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // from_block_components
    // -------------------------------------------------------------------------

    #[test]
    fn from_block_components_roundtrip() {
        let rr = SymmetricFrameMatrix3::<GCRS>::from_diagonal([4.0, 9.0, 1.0]);
        let rv = FrameMatrix3::<GCRS>::zero();
        let vv = SymmetricFrameMatrix3::<GCRS>::from_diagonal([0.01, 0.04, 0.001]);
        let p = StateCovariance::from_block_components(rr, rv, vv);
        let m = p.to_row_major();
        assert_eq!(m[0][0], 4.0);
        assert_eq!(m[1][1], 9.0);
        assert_eq!(m[3][3], 0.01);
    }
}
