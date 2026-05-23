//! # RTN orbit-to-orbit comparison
//!
//! ## Scientific scope
//!
//! Comparing an estimated trajectory against a reference orbit is most
//! interpretable in the local radial-transverse-normal frame, where along-
//! track, radial, and cross-track errors separate naturally. This module
//! computes that frame-aware difference for matched epochs.
//!
//! It assumes both trajectories are already sampled at identical epochs and
//! expressed in compatible inertial coordinates. Interpolation and frame
//! reconciliation are intentionally left to caller code.
//!
//! ## Technical scope
//!
//! The public items are `RtnDiff`, `RtnSummary`, `rtn_diff`, and
//! `rtn_summary`. Callers pass slices of `OrbitState` values and receive
//! per-epoch RTN differences plus a condensed summary for QC reporting.
//!
//! This module does not estimate the orbit or ingest external product files
//! directly.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use affn::cartesian::Displacement;
use affn::frames::GCRS;
use qtty::unit::Kilometer;
use crate::astro::dynamics::frames::RTN;
use crate::astro::dynamics::OrbitState;

/// Per-epoch RTN difference between an estimated and reference state.
#[derive(Debug, Clone)]
pub struct RtnDiff {
    /// Julian Date (TT) of the comparison epoch.
    pub jd_tt: f64,
    /// Radial difference (m).
    pub r_m: f64,
    /// Along-track difference (m).
    pub t_m: f64,
    /// Cross-track difference (m).
    pub n_m: f64,
}

impl RtnDiff {
    /// Position difference as a typed displacement in the RTN frame (kilometres).
    pub fn position_rtn_km(&self) -> Displacement<RTN, Kilometer> {
        Displacement::new(self.r_m * 1e-3, self.t_m * 1e-3, self.n_m * 1e-3)
    }
}

/// Aggregate RTN statistics over a comparison run.
#[derive(Debug, Clone)]
pub struct RtnSummary {
    /// Number of compared epochs.
    pub n: usize,
    /// 3-D RMS, metres.
    pub rms_3d_m: f64,
    /// RMS of radial component, metres.
    pub rms_r_m: f64,
    /// RMS of along-track, metres.
    pub rms_t_m: f64,
    /// RMS of cross-track, metres.
    pub rms_n_m: f64,
}

/// Pairwise RTN difference between estimated and reference states; epochs are
/// matched by index (caller is responsible for alignment and time-tag check).
pub fn rtn_diff(estimated: &[OrbitState], reference: &[OrbitState]) -> Vec<RtnDiff> {
    let n = estimated.len().min(reference.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let e = estimated[i];
        let r = reference[i];
        // Position difference in GCRS (typed Displacement<GCRS, Km>).
        let d_gcrs: Displacement<GCRS, Kilometer> = e.position - r.position;
        let frame = crate::astro::dynamics::frames::LocalOrbitalFrame::<RTN>::try_from_state(&r)
            .expect("RTN frame from reference state");
        let d_rtn = frame.to_local(d_gcrs);
        out.push(RtnDiff {
            jd_tt: e.epoch.to::<tempoch::JD>().value(),
            r_m: d_rtn.x().value() * 1000.0,
            t_m: d_rtn.y().value() * 1000.0,
            n_m: d_rtn.z().value() * 1000.0,
        });
    }
    out
}

/// Aggregate RMS statistics over a list of RTN differences.
pub fn rtn_summary(diffs: &[RtnDiff]) -> RtnSummary {
    let n = diffs.len();
    if n == 0 {
        return RtnSummary {
            n: 0,
            rms_3d_m: 0.0,
            rms_r_m: 0.0,
            rms_t_m: 0.0,
            rms_n_m: 0.0,
        };
    }
    let mut sr = 0.0;
    let mut st = 0.0;
    let mut sn = 0.0;
    for d in diffs {
        sr += d.r_m * d.r_m;
        st += d.t_m * d.t_m;
        sn += d.n_m * d.n_m;
    }
    let rms_r = (sr / n as f64).sqrt();
    let rms_t = (st / n as f64).sqrt();
    let rms_n = (sn / n as f64).sqrt();
    let rms_3d = (rms_r * rms_r + rms_t * rms_t + rms_n * rms_n).sqrt();
    RtnSummary {
        n,
        rms_3d_m: rms_3d,
        rms_r_m: rms_r,
        rms_t_m: rms_t,
        rms_n_m: rms_n,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::time::JulianDate;

    #[test]
    fn rtn_zero_for_identical_orbits() {
        let s = OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::new(7000.0, 0.0, 0.0),
            Velocity::new(0.0, 7.5, 0.0),
        );
        let d = rtn_diff(&[s], &[s]);
        assert_eq!(d.len(), 1);
        assert!(d[0].r_m.abs() < 1e-9);
        let summary = rtn_summary(&d);
        assert!(summary.rms_3d_m < 1e-9);
        assert_eq!(summary.n, 1);
    }

    #[test]
    fn radial_offset_only_appears_in_r_component() {
        let r = OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::new(7000.0, 0.0, 0.0),
            Velocity::new(0.0, 7.5, 0.0),
        );
        let e = OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::new(7000.001, 0.0, 0.0),
            Velocity::new(0.0, 7.5, 0.0),
        );
        let d = rtn_diff(&[e], &[r]);
        assert!((d[0].r_m - 1.0).abs() < 1e-9);
        assert!(d[0].t_m.abs() < 1e-9);
        assert!(d[0].n_m.abs() < 1e-9);
    }
}
