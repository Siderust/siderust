//! Declarative thrust-arc configuration that surfaces estimable parameters.
//!
//! The physical thrust model (acceleration, mass-flow rate) belongs in the
//! reusable `siderust-dynamics::thrust` module (Phase 5). This module only
//! describes a thrust arc and the parameter knobs that the estimator may
//! adjust: a scalar thrust scale and an optional 3-component delta-v at
//! arc start.

use crate::pod::problem::parameter::{Parameter, ParameterKind};

/// Declarative thrust arc.
#[derive(Debug, Clone, PartialEq)]
pub struct ThrustArcConfig {
    /// Caller-defined arc identifier.
    pub id: String,
    /// Arc start in seconds since J2000 TDB.
    pub start_seconds_tdb: f64,
    /// Arc stop in seconds since J2000 TDB.
    pub stop_seconds_tdb: f64,
    /// Estimate a multiplicative thrust-magnitude scale around 1.0.
    pub estimate_thrust_scale: bool,
    /// Estimate an unmodelled delta-v applied at arc start (m/s).
    pub estimate_delta_v: bool,
}

impl ThrustArcConfig {
    /// Emit the parameters this arc contributes to the estimator vector.
    pub fn parameters(&self) -> Vec<Parameter> {
        let mut out = Vec::new();
        if self.estimate_thrust_scale {
            out.push(Parameter {
                kind: ParameterKind::Custom(format!("thrust_scale[{}]", self.id)),
                initial_value: 1.0,
                apriori_sigma: Some(0.1),
            });
        }
        if self.estimate_delta_v {
            for axis in ['X', 'Y', 'Z'] {
                out.push(Parameter {
                    kind: ParameterKind::Custom(format!("delta_v_{axis}[{}]", self.id)),
                    initial_value: 0.0,
                    apriori_sigma: Some(0.01),
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_expected_parameters() {
        let arc = ThrustArcConfig {
            id: "burn1".into(),
            start_seconds_tdb: 0.0,
            stop_seconds_tdb: 60.0,
            estimate_thrust_scale: true,
            estimate_delta_v: true,
        };
        assert_eq!(arc.parameters().len(), 4);
    }
}
