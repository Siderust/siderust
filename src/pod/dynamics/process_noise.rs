// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Process-noise (`Q`) matrix construction for the sequential filter.
//!
//! The sequential EKF in `siderust-pod-estimation` consumes a discrete-time
//! process-noise covariance `Q(Δt)` block-diagonal in:
//!
//! 1. Cartesian position/velocity (driven by white acceleration noise).
//! 2. A drag-scale parameter `Cd_scale` modelled as a first-order
//!    Gauss–Markov process.
//! 3. An SRP-scale parameter `Crp_scale` modelled likewise.
//! 4. White empirical-acceleration coefficients (constant RTN by default).
//!
//! ## Equations
//!
//! For continuous-time white-acceleration PSD `q` (units (km/s²)²/Hz) the
//! standard van-Loan formulas yield, per axis,
//!
//! ```text
//! Q_rr = q · Δt³ / 3      Q_rv = q · Δt² / 2      Q_vv = q · Δt
//! ```
//!
//! For a Gauss–Markov scale parameter with steady-state σ and time
//! constant τ,
//!
//! ```text
//! Q_pp(Δt) = σ² · (1 − exp(−2 Δt / τ))
//! ```
//!
//! For white empirical accelerations the RTN block is `σ² · Δt` per axis.
//!
//! All time inputs are typed [`Second`] from `qtty` (the typed time
//! quantity exposed by both `qtty` and `tempoch`).

use qtty::{KmPerSecondsSquared, Second};

use super::pod_error::PodDynamicsError;

/// Continuous-time white-acceleration PSD per axis, `(km/s²)²/Hz`.
///
/// Stored as a per-axis array so radial / cross-track noise can differ when
/// the filter is rotated into a body or RTN frame downstream.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhiteAccelPsd(pub [f64; 3]);

impl WhiteAccelPsd {
    /// Build an isotropic PSD by squaring a typed `KmPerSecondsSquared`
    /// sigma (interpreted as a one-sigma white-acceleration amplitude with
    /// 1 Hz reference bandwidth — the value is squared into PSD units).
    pub fn isotropic(sigma: KmPerSecondsSquared) -> Self {
        let q = sigma.value().powi(2);
        Self([q, q, q])
    }
}

/// First-order Gauss–Markov parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussMarkovParams {
    /// Steady-state one-sigma magnitude (dimensionless for scale parameters).
    pub sigma: f64,
    /// Correlation time constant (typed seconds).
    pub tau: Second,
}

impl GaussMarkovParams {
    /// Discrete-time variance after `dt`.
    ///
    /// Returns `σ² (1 − exp(−2 Δt / τ))`.
    ///
    /// # Errors
    ///
    /// * [`PodDynamicsError::InvalidProcessNoise`] if `tau ≤ 0` or `sigma < 0`.
    pub fn variance_after(&self, dt: Second) -> Result<f64, PodDynamicsError> {
        if self.sigma < 0.0 || !self.sigma.is_finite() {
            return Err(PodDynamicsError::InvalidProcessNoise(
                "sigma must be ≥ 0 and finite",
            ));
        }
        let tau = self.tau.value();
        if tau <= 0.0 || !tau.is_finite() {
            return Err(PodDynamicsError::InvalidProcessNoise(
                "tau must be > 0 and finite",
            ));
        }
        let arg = -2.0 * dt.value() / tau;
        Ok(self.sigma.powi(2) * (1.0 - arg.exp()))
    }
}

/// Top-level process-noise specification.
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::{ProcessNoiseModel, WhiteAccelPsd, GaussMarkovParams};
/// use qtty::{KmPerSecondsSquared, Second};
/// let m = ProcessNoiseModel {
///     position_velocity: WhiteAccelPsd::isotropic(KmPerSecondsSquared::new(1e-9)),
///     drag_scale: Some(GaussMarkovParams { sigma: 0.1, tau: Second::new(3600.0) }),
///     srp_scale: None,
///     empirical_white: None,
/// };
/// let q = m.q_matrix(Second::new(60.0)).unwrap();
/// // 6-state cartesian + 1 drag-scale = 7×7
/// assert_eq!(q.len(), 7);
/// assert_eq!(q[0].len(), 7);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct ProcessNoiseModel {
    /// White-acceleration PSD driving the position/velocity 6×6 block.
    pub position_velocity: WhiteAccelPsd,
    /// Optional drag-scale Gauss–Markov term.
    pub drag_scale: Option<GaussMarkovParams>,
    /// Optional SRP-scale Gauss–Markov term.
    pub srp_scale: Option<GaussMarkovParams>,
    /// Optional white empirical-acceleration sigma (per RTN axis, treated
    /// as `σ² · Δt`). Stored as typed sigma values.
    pub empirical_white: Option<[KmPerSecondsSquared; 3]>,
}

impl ProcessNoiseModel {
    /// Total dimension `n = 6 + (drag?) + (srp?) + 3·(empirical?)`.
    pub fn dim(&self) -> usize {
        6 + usize::from(self.drag_scale.is_some())
            + usize::from(self.srp_scale.is_some())
            + if self.empirical_white.is_some() { 3 } else { 0 }
    }

    /// Construct the discrete-time process-noise covariance for the
    /// requested time step `dt`.
    ///
    /// Layout of the returned matrix is, in order:
    ///
    /// 1. `0..3` — position (km²)
    /// 2. `3..6` — velocity (km²/s²)
    /// 3. (optional) drag scale (dimensionless²)
    /// 4. (optional) SRP scale (dimensionless²)
    /// 5. (optional) 3 × empirical RTN (km/s²)²
    ///
    /// # Errors
    ///
    /// Returns [`PodDynamicsError::InvalidProcessNoise`] when any tau ≤ 0
    /// or any sigma is negative or non-finite.
    pub fn q_matrix(&self, dt: Second) -> Result<Vec<Vec<f64>>, PodDynamicsError> {
        let dt_s = dt.value();
        if dt_s < 0.0 || !dt_s.is_finite() {
            return Err(PodDynamicsError::InvalidProcessNoise(
                "dt must be ≥ 0 and finite",
            ));
        }
        let n = self.dim();
        let mut q = vec![vec![0.0_f64; n]; n];

        // Position/velocity block (white acceleration noise).
        let dt2 = dt_s * dt_s;
        let dt3 = dt2 * dt_s;
        for i in 0..3 {
            let qi = self.position_velocity.0[i];
            if qi < 0.0 || !qi.is_finite() {
                return Err(PodDynamicsError::InvalidProcessNoise(
                    "white-acceleration PSD must be ≥ 0 and finite",
                ));
            }
            q[i][i] = qi * dt3 / 3.0;
            q[3 + i][3 + i] = qi * dt_s;
            let cross = qi * dt2 / 2.0;
            q[i][3 + i] = cross;
            q[3 + i][i] = cross;
        }

        let mut idx = 6;
        if let Some(gm) = self.drag_scale {
            q[idx][idx] = gm.variance_after(dt)?;
            idx += 1;
        }
        if let Some(gm) = self.srp_scale {
            q[idx][idx] = gm.variance_after(dt)?;
            idx += 1;
        }
        if let Some(emp) = self.empirical_white {
            for k in 0..3 {
                let s = emp[k].value();
                if s < 0.0 || !s.is_finite() {
                    return Err(PodDynamicsError::InvalidProcessNoise(
                        "empirical sigma must be ≥ 0 and finite",
                    ));
                }
                q[idx + k][idx + k] = s * s * dt_s;
            }
        }
        Ok(q)
    }

    /// Return `true` iff the matrix produced by `q_matrix(dt)` is positive
    /// semi-definite. Implemented by Cholesky-with-pivot of the symmetric
    /// matrix; for the small dimensions we deal with (≤ 12) it is exact and
    /// allocation-free w.r.t. the result, and avoids pulling in a heavy
    /// dependency.
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::{ProcessNoiseModel, WhiteAccelPsd, GaussMarkovParams};
    /// use qtty::{KmPerSecondsSquared, Second};
    /// let m = ProcessNoiseModel {
    ///     position_velocity: WhiteAccelPsd::isotropic(KmPerSecondsSquared::new(1e-9)),
    ///     drag_scale: Some(GaussMarkovParams { sigma: 0.1, tau: Second::new(3600.0) }),
    ///     srp_scale: None,
    ///     empirical_white: None,
    /// };
    /// assert!(m.is_psd(Second::new(60.0)).unwrap());
    /// ```
    pub fn is_psd(&self, dt: Second) -> Result<bool, PodDynamicsError> {
        let q = self.q_matrix(dt)?;
        Ok(matrix_is_psd(&q))
    }
}

/// Test a small symmetric matrix for positive semi-definiteness via a
/// pivoted Cholesky-style factorisation with a small absolute tolerance.
fn matrix_is_psd(a: &[Vec<f64>]) -> bool {
    let n = a.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    let tol = 1e-18;
    for i in 0..n {
        let mut s = a[i][i];
        s -= l[i][..i].iter().map(|x| x * x).sum::<f64>();
        if s < -tol {
            return false;
        }
        let lii = s.max(0.0).sqrt();
        l[i][i] = lii;
        if lii == 0.0 {
            // Whole row must be zero for PSD with zero pivot.
            for j in (i + 1)..n {
                let s2 = a[j][i]
                    - l[j][..i]
                        .iter()
                        .zip(l[i][..i].iter())
                        .map(|(a, b)| a * b)
                        .sum::<f64>();
                if s2.abs() > 1e-12 {
                    return false;
                }
                l[j][i] = 0.0;
            }
        } else {
            for j in (i + 1)..n {
                let s2 = a[j][i]
                    - l[j][..i]
                        .iter()
                        .zip(l[i][..i].iter())
                        .map(|(a, b)| a * b)
                        .sum::<f64>();
                l[j][i] = s2 / lii;
            }
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessNoise enum (POD-layer high-level API)
// ─────────────────────────────────────────────────────────────────────────────

/// One interval in a piecewise-constant noise schedule.
///
/// Units are **SI** (meters / m/s) rather than km / km/s so that the values
/// align with typical OD covariance budgets expressed in SI.
///
/// # Example
///
/// ```
/// use qtty::Second;
/// use siderust::pod::dynamics::process_noise::PiecewiseSegment;
///
/// let seg = PiecewiseSegment {
///     duration: Second::new(300.0),
///     sigma_pos_m: 1.0,
///     sigma_vel_mps: 0.001,
/// };
/// assert_eq!(seg.duration.value(), 300.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PiecewiseSegment {
    /// Duration of this noise interval (seconds).
    pub duration: Second,
    /// One-sigma position noise amplitude (metres).
    pub sigma_pos_m: f64,
    /// One-sigma velocity noise amplitude (metres per second).
    pub sigma_vel_mps: f64,
}

/// High-level process-noise specification for a single filter update step.
///
/// Wraps the common noise configurations in a convenient enum so that
/// downstream filter code can dispatch without inspecting struct fields.
/// Differs from [`ProcessNoiseModel`] in that:
///
/// * Units are SI (metres, m/s) rather than km, km/s.
/// * The type is an `enum` for ergonomic matching rather than a struct.
/// * Drag/SRP scale and empirical-acceleration terms are handled by
///   [`ProcessNoiseModel`]; this enum covers the position/velocity noise budget.
///
/// # Example
///
/// ```
/// use siderust::pod::dynamics::process_noise::{PiecewiseSegment, ProcessNoise};
/// use qtty::Second;
///
/// // No process noise at all.
/// let n = ProcessNoise::None;
/// assert!(!n.is_active());
///
/// // Isotropic white noise.
/// let w = ProcessNoise::WhiteNoise { sigma_pos_m: 10.0, sigma_vel_mps: 0.01 };
/// assert!(w.is_active());
///
/// // Two-segment piecewise-constant schedule.
/// let p = ProcessNoise::PiecewiseConstant {
///     intervals: vec![
///         PiecewiseSegment { duration: Second::new(60.0), sigma_pos_m: 5.0, sigma_vel_mps: 0.005 },
///         PiecewiseSegment { duration: Second::new(300.0), sigma_pos_m: 10.0, sigma_vel_mps: 0.01 },
///     ],
/// };
/// assert!(p.is_active());
/// if let ProcessNoise::PiecewiseConstant { intervals } = &p {
///     assert_eq!(intervals.len(), 2);
/// }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum ProcessNoise {
    /// No process noise: the covariance is evolved by the STM only.
    None,

    /// Isotropic white-noise injection with fixed sigma values for all steps.
    ///
    /// Units: metres (position) and metres per second (velocity).
    WhiteNoise {
        /// One-sigma position noise (metres).
        sigma_pos_m: f64,
        /// One-sigma velocity noise (metres per second).
        sigma_vel_mps: f64,
    },

    /// Piecewise-constant noise schedule: a different sigma level applies to
    /// each time interval. The interval at index `i` applies for a step whose
    /// elapsed arc time falls within `sum(durations[0..i])..sum(durations[0..=i])`.
    ///
    /// Intervals are matched by cumulative elapsed arc time. If the current
    /// time lies beyond the last interval, the last segment's sigma values are
    /// used (constant extension).
    PiecewiseConstant {
        /// Ordered list of intervals (first interval applies earliest).
        intervals: Vec<PiecewiseSegment>,
    },
}

impl ProcessNoise {
    /// `true` iff this variant injects any noise (i.e. is not [`ProcessNoise::None`]).
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::process_noise::ProcessNoise;
    /// assert!(!ProcessNoise::None.is_active());
    /// assert!(ProcessNoise::WhiteNoise { sigma_pos_m: 1.0, sigma_vel_mps: 0.01 }.is_active());
    /// ```
    pub fn is_active(&self) -> bool {
        !matches!(self, ProcessNoise::None)
    }

    /// Look up the noise sigmas applicable at cumulative elapsed time `elapsed`.
    ///
    /// Returns `(sigma_pos_m, sigma_vel_mps)`:
    ///
    /// * [`ProcessNoise::None`] → `(0.0, 0.0)`.
    /// * [`ProcessNoise::WhiteNoise`] → the fixed sigma pair.
    /// * [`ProcessNoise::PiecewiseConstant`] → sigma of the interval that
    ///   contains `elapsed`; falls back to the last interval for times beyond
    ///   the schedule. Returns `(0.0, 0.0)` for an empty schedule.
    ///
    /// # Example
    ///
    /// ```
    /// use siderust::pod::dynamics::process_noise::{PiecewiseSegment, ProcessNoise};
    /// use qtty::Second;
    ///
    /// let p = ProcessNoise::PiecewiseConstant {
    ///     intervals: vec![
    ///         PiecewiseSegment { duration: Second::new(60.0), sigma_pos_m: 5.0, sigma_vel_mps: 0.005 },
    ///         PiecewiseSegment { duration: Second::new(300.0), sigma_pos_m: 10.0, sigma_vel_mps: 0.01 },
    ///     ],
    /// };
    /// // Before first boundary (30 s < 60 s).
    /// assert_eq!(p.sigma_at(Second::new(30.0)), (5.0, 0.005));
    /// // After first boundary (100 s in 60..360 range).
    /// assert_eq!(p.sigma_at(Second::new(100.0)), (10.0, 0.01));
    /// // Beyond schedule → last segment.
    /// assert_eq!(p.sigma_at(Second::new(1000.0)), (10.0, 0.01));
    /// ```
    pub fn sigma_at(&self, elapsed: Second) -> (f64, f64) {
        match self {
            ProcessNoise::None => (0.0, 0.0),
            ProcessNoise::WhiteNoise {
                sigma_pos_m,
                sigma_vel_mps,
            } => (*sigma_pos_m, *sigma_vel_mps),
            ProcessNoise::PiecewiseConstant { intervals } => {
                if intervals.is_empty() {
                    return (0.0, 0.0);
                }
                let mut cumulative = 0.0;
                for seg in intervals {
                    cumulative += seg.duration.value();
                    if elapsed.value() < cumulative {
                        return (seg.sigma_pos_m, seg.sigma_vel_mps);
                    }
                }
                // Elapsed is past the last boundary — use last segment.
                let last = intervals.last().expect("non-empty checked above");
                (last.sigma_pos_m, last.sigma_vel_mps)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimensions_track_options() {
        let m = ProcessNoiseModel {
            position_velocity: WhiteAccelPsd::isotropic(KmPerSecondsSquared::new(1e-9)),
            drag_scale: Some(GaussMarkovParams {
                sigma: 0.1,
                tau: Second::new(3600.0),
            }),
            srp_scale: Some(GaussMarkovParams {
                sigma: 0.05,
                tau: Second::new(7200.0),
            }),
            empirical_white: Some([KmPerSecondsSquared::new(1e-10); 3]),
        };
        assert_eq!(m.dim(), 6 + 1 + 1 + 3);
        let q = m.q_matrix(Second::new(60.0)).unwrap();
        assert_eq!(q.len(), 11);
    }

    #[test]
    fn rejects_negative_sigma() {
        let m = ProcessNoiseModel {
            position_velocity: WhiteAccelPsd::isotropic(KmPerSecondsSquared::new(1e-9)),
            drag_scale: Some(GaussMarkovParams {
                sigma: -1.0,
                tau: Second::new(60.0),
            }),
            srp_scale: None,
            empirical_white: None,
        };
        assert!(matches!(
            m.q_matrix(Second::new(10.0)).unwrap_err(),
            PodDynamicsError::InvalidProcessNoise(_)
        ));
    }

    #[test]
    fn process_noise_none_is_inactive() {
        assert!(!ProcessNoise::None.is_active());
        assert_eq!(ProcessNoise::None.sigma_at(Second::new(0.0)), (0.0, 0.0));
    }

    #[test]
    fn process_noise_white_noise_is_active() {
        let n = ProcessNoise::WhiteNoise {
            sigma_pos_m: 5.0,
            sigma_vel_mps: 0.01,
        };
        assert!(n.is_active());
        assert_eq!(n.sigma_at(Second::new(100.0)), (5.0, 0.01));
    }

    #[test]
    fn piecewise_constant_dispatches_by_elapsed_time() {
        let p = ProcessNoise::PiecewiseConstant {
            intervals: vec![
                PiecewiseSegment {
                    duration: Second::new(60.0),
                    sigma_pos_m: 1.0,
                    sigma_vel_mps: 0.001,
                },
                PiecewiseSegment {
                    duration: Second::new(300.0),
                    sigma_pos_m: 10.0,
                    sigma_vel_mps: 0.01,
                },
            ],
        };
        assert_eq!(p.sigma_at(Second::new(0.0)), (1.0, 0.001));
        assert_eq!(p.sigma_at(Second::new(59.9)), (1.0, 0.001));
        assert_eq!(p.sigma_at(Second::new(60.0)), (10.0, 0.01));
        // Beyond schedule → last segment.
        assert_eq!(p.sigma_at(Second::new(9_999.0)), (10.0, 0.01));
    }

    #[test]
    fn piecewise_constant_empty_schedule_returns_zero() {
        let p = ProcessNoise::PiecewiseConstant { intervals: vec![] };
        assert_eq!(p.sigma_at(Second::new(0.0)), (0.0, 0.0));
    }
}
