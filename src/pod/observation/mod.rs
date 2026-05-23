//! # siderust-pod observation models
//!
//! ## Scientific scope
//!
//! This module houses the measurement equations that connect estimated
//! spacecraft states to tracked observables. The supported observation
//! types are:
//!
//! - **GNSS** pseudorange and carrier-phase (Saastamoinen troposphere,
//!   Klobuchar ionosphere, satellite/receiver clock providers, Sagnac
//!   and Shapiro corrections).
//! - **SLR** two-way normal-point range (Marini-Murray troposphere,
//!   Shapiro relativistic delay, range-bias support).
//! - Generic typed [`Observation`] trait for
//!   object-safe batching plus a [`CorrectionRegistry`] pipeline
//!   (PCO/PCV, Shapiro, Earth tide displacement).
//!
//! Mission-specific observation models (e.g. LISA inter-satellite range)
//! are implemented in `siderust/examples/` rather than the library source.
//!
//! ## Technical scope
//!
//! Every observation type is parameterized over its provider bundle
//! (clock, atmosphere, antenna), accepts typed states
//! (`Time<S, F>` epochs, `Position` and `Velocity` in `affn` frames),
//! and returns scalar predictions plus partial derivatives ready for
//! the estimation layer. Parsing of raw data files and orchestration of
//! estimation loops are out of scope here.
//!
//! ## References
//!
//! - Misra, P., & Enge, P. (2012). *Global Positioning System: Signals,
//!   Measurements, and Performance* (2nd ed.). Ganga-Jamuna Press.
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical
//!   Orbit Determination*. Elsevier Academic Press.
//! - International Laser Ranging Service (2009). *ILRS Standard Data
//!   Formats*.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod batch;
pub mod corrections;
pub mod error;
pub mod gnss_obs;
pub mod obs_trait;
pub mod provider_bundle;
pub mod slr_obs;

pub use batch::ObservationBatch;
pub use corrections::{
    Correction, CorrectionRegistry, EarthTideDisplacement, PhaseCenterOffset, ShapiroDelay,
};
pub use error::PodObservationsError;
pub use gnss_obs::{
    GnssCarrierPhaseObs, GnssPseudorangeObs, IonoModel, KlobucharParams, TropModel,
};
pub use obs_trait::{
    AnyObservation, CartesianState, ObsResidual, ObsType, Observation, PhaseResidual,
};
pub use provider_bundle::{NullProviderBundle, ProviderBundle};
pub use slr_obs::SlrNormalPointObs;

