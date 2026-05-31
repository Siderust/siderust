//! Parameter typing for the POD estimator.
//!
//! The estimator solves for a heterogeneous parameter vector: orbital
//! state components, receiver clock, drag scale, SRP scale, carrier float
//! ambiguities, range biases, etc. [`ParameterKind`] tags each scalar slot
//! so that downstream code (covariance ordering, residual reports, manifest
//! dumps) preserves semantic meaning.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Semantic kind of a single estimated parameter.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ParameterKind {
    /// X component of Cartesian position at a reference epoch.
    StatePositionX,
    /// Y component of Cartesian position at a reference epoch.
    StatePositionY,
    /// Z component of Cartesian position at a reference epoch.
    StatePositionZ,
    /// X component of Cartesian velocity at a reference epoch.
    StateVelocityX,
    /// Y component of Cartesian velocity at a reference epoch.
    StateVelocityY,
    /// Z component of Cartesian velocity at a reference epoch.
    StateVelocityZ,
    /// Receiver clock bias (seconds) for a given station/receiver id.
    ReceiverClockBias {
        /// Receiver or station identifier.
        receiver_id: String,
    },
    /// Drag scale factor (dimensionless multiplier on the nominal drag force).
    DragScale,
    /// Solar-radiation-pressure scale factor.
    SrpScale,
    /// Carrier-phase float ambiguity for a (receiver, satellite, signal) triple, in cycles.
    CarrierAmbiguity {
        /// Receiver or station identifier.
        receiver_id: String,
        /// Satellite PRN or identifier.
        satellite_id: String,
        /// Signal band/code designation (e.g., `"L1"`).
        signal: String,
    },
    /// Range bias for an SLR station, in meters.
    SlrRangeBias {
        /// SLR station identifier.
        station_id: String,
    },
    /// Empirical acceleration in a local frame (e.g. RTN) component.
    EmpiricalAcceleration {
        /// Frame identifier (e.g., `"RTN"`).
        frame: String,
        /// Component axis (`'R'`, `'T'`, `'N'`, etc.).
        component: char,
    },
    /// Custom user-defined parameter for extension models.
    Custom(String),
}

/// One slot in the estimator's parameter vector.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Parameter {
    /// Semantic kind of this parameter.
    pub kind: ParameterKind,
    /// Initial value used by the estimator (engineering units appropriate to `kind`).
    pub initial_value: f64,
    /// Optional a-priori 1-sigma constraint; `None` means unconstrained.
    pub apriori_sigma: Option<f64>,
}

/// Total parameter ordering for a run; the estimator must use this to lay
/// out residual rows and design columns.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ParameterOrdering {
    /// List of parameters in the order used by the estimator.
    pub params: Vec<Parameter>,
}

impl ParameterOrdering {
    /// Number of parameters in the ordering.
    pub fn len(&self) -> usize {
        self.params.len()
    }
    /// Whether the parameter list is empty.
    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }
    /// Index of the first parameter matching `kind`, if any.
    pub fn index_of(&self, kind: &ParameterKind) -> Option<usize> {
        self.params.iter().position(|p| &p.kind == kind)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordering_lookup() {
        let ord = ParameterOrdering {
            params: vec![
                Parameter {
                    kind: ParameterKind::StatePositionX,
                    initial_value: 7000.0,
                    apriori_sigma: Some(0.1),
                },
                Parameter {
                    kind: ParameterKind::DragScale,
                    initial_value: 1.0,
                    apriori_sigma: Some(0.2),
                },
            ],
        };
        assert_eq!(ord.index_of(&ParameterKind::DragScale), Some(1));
        assert_eq!(ord.index_of(&ParameterKind::SrpScale), None);
    }
}
