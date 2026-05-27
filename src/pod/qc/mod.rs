//! # siderust::pod quality control
//!
//! ## Scientific scope
//!
//! Post-fit diagnostics for POD workflows: residual summaries,
//! orbit-to-orbit comparisons in local frames, and SLR validation summaries.
//! The intent is to quantify solution quality without modifying the estimated
//! orbit itself.
//!
//! ## Technical scope
//!
//! Orbit-comparison helpers, grouped residual statistics, and SLR validation
//! report types. Inputs are already-computed trajectories, residual lists, or
//! QC JSON payloads prepared by upstream stages.
//!
//! HTML rendering of QC reports is a service concern and lives in
//! `constops::report::html`.
//!
//! This module does not read raw measurement files or perform estimation.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
#![forbid(unsafe_code)]

pub mod orbit_compare;
pub mod residuals;
pub mod slr_validation;

pub use orbit_compare::{rtn_diff, rtn_summary, RtnDiff, RtnSummary};
pub use residuals::{ResidualStats, ResidualsByGroup};
pub use slr_validation::{SlrResidual, SlrValidationReport};
