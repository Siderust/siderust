//! # Residual statistics
//!
//! ## Scientific scope
//!
//! Residual statistics provide the first compact view of whether a POD
//! solution is fitting its measurements consistently. This module groups
//! scalar residuals and derives summary quantities such as RMS, mean, and
//! extrema for QC reporting.
//!
//! The statistics are classical descriptive measures; there is no outlier
//! rejection or stochastic modelling in this layer. Their interpretation
//! depends on the observation model and data editing done upstream.
//!
//! ## Technical scope
//!
//! The public types are `ResidualStats` and `ResidualsByGroup`, which
//! organize scalar residual collections into report-friendly summaries.
//! Callers feed already-computed residual values and group labels into this
//! layer.
//!
//! Residual formation, parameter estimation, and product serialization are
//! outside its scope.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Summary statistics for a single residual stream.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ResidualStats {
    /// Number of residuals.
    pub n: usize,
    /// Mean (signed) residual.
    pub mean: f64,
    /// Standard deviation (sample, n-1 denominator if n > 1).
    pub std: f64,
    /// Root-mean-square value.
    pub rms: f64,
    /// Min value.
    pub min: f64,
    /// Max value.
    pub max: f64,
}

impl ResidualStats {
    /// Compute statistics for a slice.
    pub fn from_slice(values: &[f64]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::default();
        }
        let sum: f64 = values.iter().sum();
        let mean = sum / n as f64;
        let mut var = 0.0;
        let mut sq = 0.0;
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for &v in values {
            let d = v - mean;
            var += d * d;
            sq += v * v;
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }
        let std = if n > 1 {
            (var / (n as f64 - 1.0)).sqrt()
        } else {
            0.0
        };
        let rms = (sq / n as f64).sqrt();
        Self {
            n,
            mean,
            std,
            rms,
            min,
            max,
        }
    }
}

/// Stats grouped by an arbitrary string key (e.g. measurement type).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ResidualsByGroup {
    /// Per-group statistics.
    pub groups: BTreeMap<String, ResidualStats>,
    /// Overall stats across every group.
    pub overall: ResidualStats,
}

impl ResidualsByGroup {
    /// Compute grouped + overall stats from `(group, value)` pairs.
    pub fn from_pairs<'a, I>(items: I) -> Self
    where
        I: IntoIterator<Item = (&'a str, f64)>,
    {
        let mut buckets: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        let mut all = Vec::new();
        for (g, v) in items {
            buckets.entry(g.to_string()).or_default().push(v);
            all.push(v);
        }
        let groups = buckets
            .into_iter()
            .map(|(k, v)| (k, ResidualStats::from_slice(&v)))
            .collect();
        Self {
            groups,
            overall: ResidualStats::from_slice(&all),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_stats() {
        let s = ResidualStats::from_slice(&[1.0, -1.0, 1.0, -1.0]);
        assert_eq!(s.n, 4);
        assert!((s.mean).abs() < 1e-12);
        assert!((s.rms - 1.0).abs() < 1e-12);
    }

    #[test]
    fn grouped() {
        let g = ResidualsByGroup::from_pairs(vec![
            ("code", 0.5),
            ("code", -0.5),
            ("phase", 0.01),
            ("phase", -0.01),
        ]);
        assert_eq!(g.groups.len(), 2);
        assert!((g.groups["code"].rms - 0.5).abs() < 1e-12);
        assert!((g.overall.n) == 4);
    }
}
