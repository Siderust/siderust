//! # Sequential Kalman estimation
//!
//! ## Scientific scope
//!
//! This module provides a sequential estimator shaped like an Extended
//! Kalman Filter for POD-style state updates. It is intended for settings
//! where observations arrive one at a time and the caller can supply both a
//! transition model and scalar measurement linearization.
//!
//! The implementation targets compact, deterministic replay scenarios.
//! Process-noise modelling, RTS smoothing, and advanced numerical
//! stabilization are intentionally kept minimal at this stage and are
//! tracked in the long-term roadmap (`reqs-and-plan.md`).
//!
//! ## Technical scope
//!
//! The main public items are `OrbitEkf`, `EkfError`, and `InnovationRecord`.
//! `OrbitEkf` is specialised for 6D orbit POD and stores [`OrbitState`] and
//! [`StateCovariance<GCRS>`] as first-class fields.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use crate::astro::dynamics::{OrbitState, Velocity};
use crate::coordinates::frames::GCRS;
use affn::Displacement;
use faer::Mat;
use principia::StateCovariance;
use qtty::dynamics::KmPerSeconds;
use qtty::length::Kilometers;
use qtty::unit::Kilometer;
use thiserror::Error;

/// EKF error type.
#[derive(Debug, Error)]
pub enum EkfError {
    /// Singular innovation covariance (S ≤ 0 for a scalar measurement).
    #[error("singular innovation covariance: {0}")]
    Singular(f64),
}

/// Per-measurement innovation/NIS record produced by `update_scalar`.
#[derive(Debug, Clone, Copy)]
pub struct InnovationRecord {
    /// Measurement minus prediction.
    pub innovation: f64,
    /// Innovation variance S = HPHᵀ + R.
    pub variance: f64,
    /// Normalised innovation squared (innovation² / S).
    pub nis: f64,
}

/// Typed 6D orbit EKF specialised for position + velocity state.
///
/// Stores the filter state directly as siderust/affn types:
///
/// - [`OrbitState`] for the mean (position, velocity, epoch)
/// - [`StateCovariance<GCRS>`] for the 6×6 covariance split into typed 3×3
///   blocks
///
/// All internal linear algebra is performed on a `faer::Mat<f64>` extracted
/// from the typed containers via [`StateCovariance::to_row_major`] and packed
/// back with [`StateCovariance::from_blocks`].
#[derive(Debug, Clone)]
pub struct OrbitEkf {
    state: OrbitState,
    cov: StateCovariance<GCRS>,
}

impl OrbitEkf {
    /// New filter with explicit typed state and covariance.
    pub fn new(state: OrbitState, cov: StateCovariance<GCRS>) -> Self {
        Self { state, cov }
    }

    /// New filter with diagonal covariance built from per-axis 1-σ position
    /// (km) and velocity (km/s) standard deviations.
    pub fn from_stddevs(state: OrbitState, sigma_pos: [f64; 3], sigma_vel: [f64; 3]) -> Self {
        let sigma_pos_typed = [
            Kilometers::new(sigma_pos[0]),
            Kilometers::new(sigma_pos[1]),
            Kilometers::new(sigma_pos[2]),
        ];
        let sigma_vel_typed = [
            KmPerSeconds::new(sigma_vel[0]),
            KmPerSeconds::new(sigma_vel[1]),
            KmPerSeconds::new(sigma_vel[2]),
        ];
        Self {
            state,
            cov: StateCovariance::<GCRS>::diagonal_from_sigmas(sigma_pos_typed, sigma_vel_typed),
        }
    }

    /// Current orbit state (position + velocity + epoch).
    pub fn state(&self) -> &OrbitState {
        &self.state
    }

    /// Current 6×6 covariance in GCRS.
    pub fn covariance(&self) -> &StateCovariance<GCRS> {
        &self.cov
    }

    /// Time update: replace the mean with `state_pred` and propagate covariance
    /// as `P ← Φ P Φᵀ + Q`.
    ///
    /// `phi` is the 6×6 state-transition matrix in row-major order.
    /// `q` is an optional process-noise covariance added after the propagation.
    pub fn predict(
        &mut self,
        state_pred: OrbitState,
        phi: [[f64; 6]; 6],
        q: Option<StateCovariance<GCRS>>,
    ) {
        let phi_m: Mat<f64> = Mat::from_fn(6, 6, |i, j| phi[i][j]);
        let p = self.cov.to_row_major();
        let p_m: Mat<f64> = Mat::from_fn(6, 6, |i, j| p[i][j]);
        // P ← Φ P Φᵀ (+ Q)
        let phi_p = &phi_m * &p_m;
        let mut p_next = &phi_p * phi_m.transpose();
        if let Some(q_cov) = q {
            let q = q_cov.to_row_major();
            let q_m: Mat<f64> = Mat::from_fn(6, 6, |i, j| q[i][j]);
            p_next += q_m;
        }
        self.state = state_pred;
        self.cov = StateCovariance::from_row_major(std::array::from_fn(|i| {
            std::array::from_fn(|j| p_next[(i, j)])
        }));
    }

    /// Scalar measurement update.
    ///
    /// `h` is the 6-element row of partial derivatives ∂y/∂x.
    /// `r` is the measurement variance σ².
    /// Returns the innovation record or an error if the innovation variance
    /// is not strictly positive.
    pub fn update_scalar(
        &mut self,
        h: [f64; 6],
        innovation: f64,
        r: f64,
    ) -> Result<InnovationRecord, EkfError> {
        let p = self.cov.to_row_major();
        let p_m: Mat<f64> = Mat::from_fn(6, 6, |i, j| p[i][j]);
        let h_col: Mat<f64> = Mat::from_fn(6, 1, |i, _| h[i]);

        // P h
        let ph = &p_m * &h_col;

        // S = hᵀ P h + R
        let mut s = r;
        for i in 0..6 {
            s += h[i] * ph[(i, 0)];
        }
        if s.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return Err(EkfError::Singular(s));
        }

        // K = P h / S
        let k: Mat<f64> = Mat::from_fn(6, 1, |i, _| ph[(i, 0)] / s);

        // x ← x + K * innovation  (typed arithmetic)
        let pos_delta = Displacement::<GCRS, Kilometer>::new(
            k[(0, 0)] * innovation,
            k[(1, 0)] * innovation,
            k[(2, 0)] * innovation,
        );
        let vel_delta = Velocity::<GCRS>::new(
            k[(3, 0)] * innovation,
            k[(4, 0)] * innovation,
            k[(5, 0)] * innovation,
        );
        let new_pos = self.state.position + pos_delta;
        let new_vel = self.state.velocity + vel_delta;
        self.state = OrbitState::new(self.state.epoch, new_pos, new_vel);

        // P ← P − K (hᵀ P)  then symmetrise
        let htp: Mat<f64> =
            Mat::from_fn(1, 6, |_, j| (0..6).map(|i| h[i] * p_m[(i, j)]).sum::<f64>());
        let mut p_next = p_m - &k * &htp;
        // Symmetrise
        for i in 0..6 {
            for j in (i + 1)..6 {
                let v = 0.5 * (p_next[(i, j)] + p_next[(j, i)]);
                p_next[(i, j)] = v;
                p_next[(j, i)] = v;
            }
        }
        self.cov = StateCovariance::from_row_major(std::array::from_fn(|i| {
            std::array::from_fn(|j| p_next[(i, j)])
        }));
        let nis = innovation * innovation / s;
        Ok(InnovationRecord {
            innovation,
            variance: s,
            nis,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::Position;
    use crate::time::JulianDate;

    fn make_orbit_state() -> OrbitState {
        let epoch = JulianDate::new(2_451_545.0);
        let pos = Position::<GCRS>::new(7000.0, 0.0, 0.0);
        let vel = Velocity::<GCRS>::new(0.0, 7.5, 0.0);
        OrbitState::new(epoch.to_j2000s(), pos, vel)
    }

    #[test]
    fn orbit_ekf_update_reduces_variance() {
        let s0 = make_orbit_state();
        let mut f = OrbitEkf::from_stddevs(s0, [1.0, 1.0, 1.0], [1e-3, 1e-3, 1e-3]);
        let p_before = f.covariance().to_row_major();
        // Observe x-component directly (h = [1,0,0,0,0,0]).
        let rec = f
            .update_scalar([1.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.5, 0.01)
            .unwrap();
        let p_after = f.covariance().to_row_major();
        assert!(p_after[0][0] < p_before[0][0], "pos variance must shrink");
        assert!(rec.nis > 0.0);
        // State should move toward measurement.
        assert!(f.state().position.x().value() > 7000.0);
    }

    #[test]
    fn orbit_ekf_predict_propagates_covariance() {
        let s0 = make_orbit_state();
        let mut f = OrbitEkf::from_stddevs(s0, [0.1, 0.1, 0.1], [1e-4, 1e-4, 1e-4]);
        let p_before = f.covariance().to_row_major();
        // Identity STM: covariance stays the same without Q.
        let mut phi = [[0.0f64; 6]; 6];
        for i in 0..6 {
            phi[i][i] = 1.0;
        }
        f.predict(s0, phi, None);
        let p_after_no_q = f.covariance().to_row_major();
        for i in 0..6 {
            assert!(
                (p_after_no_q[i][i] - p_before[i][i]).abs() < 1e-12,
                "identity phi should leave diagonal unchanged"
            );
        }
        // Add process noise: diagonal must inflate.
        let q_pos = [
            Kilometers::new(0.01),
            Kilometers::new(0.01),
            Kilometers::new(0.01),
        ];
        let q_vel = [
            KmPerSeconds::new(1e-5),
            KmPerSeconds::new(1e-5),
            KmPerSeconds::new(1e-5),
        ];
        let q = StateCovariance::<GCRS>::diagonal_from_sigmas(q_pos, q_vel);
        f.predict(s0, phi, Some(q));
        let p_inflated = f.covariance().to_row_major();
        assert!(
            p_inflated[0][0] > p_after_no_q[0][0],
            "Q must inflate pos variance"
        );
        assert!(
            p_inflated[3][3] > p_after_no_q[3][3],
            "Q must inflate vel variance"
        );
    }
}
