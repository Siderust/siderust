//! # SLR validation summaries
//!
//! ## Scientific scope
//!
//! SLR serves as an external geometric check on an orbit solution because
//! it observes two-way range using instrumentation independent from GNSS
//! tracking. This module packages SLR observed-minus-computed residuals
//! into a compact validation report for the current workspace.
//!
//! Its regime is intentionally lightweight and follows the simplified SLR
//! measurement model implemented upstream. Operational-grade bias
//! calibration, atmospheric correction modelling, and station eccentricity
//! handling are beyond the present scope.
//!
//! ## Technical scope
//!
//! The public types are `SlrResidual` and `SlrValidationReport`. Callers
//! provide already modelled residual values and receive a report object
//! ready for QC JSON output or test assertions.
//!
//! Numerical light-time modelling is delegated to
//! `siderust::pod::observation::slr_obs`.
//!
//! ## References
//!
//! - Degnan, J. J. (1993). Millimeter accuracy satellite laser ranging: a
//!   review. Contributions of Space Geodesy to Geodynamics: Technology, 25,
//!   133-162.
//! - Pearlman, M. R., Noll, C. E., et al. (2019). The ILRS: Current status
//!   and future prospects. Journal of Geodesy, 93, 2161-2180.
use super::residuals::ResidualStats;

/// One O−C residual at an SLR observation epoch.
#[derive(Debug, Clone)]
pub struct SlrResidual {
    /// Julian Date (TT) of the bounce epoch.
    pub jd_tt: f64,
    /// O−C residual, metres.
    pub residual_m: f64,
}

/// Aggregate SLR validation report.
#[derive(Debug, Clone)]
pub struct SlrValidationReport {
    /// Per-pass residuals.
    pub residuals: Vec<SlrResidual>,
    /// Aggregate statistics.
    pub stats: ResidualStats,
}

impl SlrValidationReport {
    /// Build the report from a list of `(jd_tt, residual_m)` pairs.
    pub fn from_pairs<I: IntoIterator<Item = (f64, f64)>>(pairs: I) -> Self {
        let residuals: Vec<SlrResidual> = pairs
            .into_iter()
            .map(|(t, r)| SlrResidual {
                jd_tt: t,
                residual_m: r,
            })
            .collect();
        let values: Vec<f64> = residuals.iter().map(|r| r.residual_m).collect();
        let stats = ResidualStats::from_slice(&values);
        Self { residuals, stats }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aggregates_zero_residuals() {
        let r = SlrValidationReport::from_pairs((0..5).map(|i| (i as f64, 0.0)));
        assert_eq!(r.residuals.len(), 5);
        assert!(r.stats.rms.abs() < 1e-12);
    }
}
