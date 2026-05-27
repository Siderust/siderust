#![allow(clippy::needless_range_loop, clippy::inconsistent_digit_grouping)]
//! # siderust::pod estimation
//!
//! ## Scientific scope
//!
//! This crate collects the numerical estimation kernels used by the POD
//! workspace: batch weighted least squares, nonlinear Gauss-Newton
//! iteration, and a sequential Kalman-style update path. The scientific
//! target is short-arc orbit determination where measurement models and
//! force models are supplied by sibling modules.
//!
//! These routines assume callers provide already-linearized residual rows
//! or propagation closures at the appropriate orbit-determination epoch.
//! They do not introduce their own measurement corrections, force laws, or
//! time-scale conversions.
//!
//! ## Technical scope
//!
//! The module re-exports `gauss_newton`, `NormalEquations`, `OrbitEkf`, and the
//! parameter descriptors that define estimator-state ordering. Public APIs
//! mainly consume scalar residuals, Jacobians, and covariance-like matrices
//! in solver space, while typed orbital and temporal quantities remain
//! owned by `siderust` and the observation/service crates.
//!
//! This module is the numerical core only: file formats, synthetic data
//! generation, and product writing are intentionally delegated to sibling
//! crates.
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). Statistical Orbit
//!   Determination. Elsevier Academic Press.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
#![forbid(unsafe_code)]

pub mod nonlinear;
pub mod sequential;
pub mod wls;

pub use crate::pod::problem::parameter::{Parameter, ParameterKind};
pub use nonlinear::{gauss_newton, NonlinearError, NonlinearOptions, NonlinearReport};
pub use sequential::{EkfError, InnovationRecord, OrbitEkf};
pub use wls::{NormalEquations, WlsResult, WlsSolverError};
