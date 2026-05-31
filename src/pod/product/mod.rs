//! # siderust::pod product writers
//!
//! ## Scientific scope
//!
//! This crate turns in-memory estimation results into interchange and
//! reporting artifacts. The scientific content is inherited from the
//! upstream orbit solution and residual statistics; this layer focuses on
//! packaging those results into standard or workspace-defined product
//! shapes.
//!
//! Current outputs cover precise orbit histories, CCSDS ephemerides,
//! grouped QC JSON, tabular residual exports, Parquet residual tables
//! (feature-gated), and run manifests.
//!
//! ## Technical scope
//!
//! | Module | Surface |
//! |---|---|
//! | [`orbit`] | `write_sp3_from_states`, `write_oem_from_states`, `write_oem_from_spacecraft_states`, `Sp3ProductWriter`, `OemProductWriter` |
//! | [`naming`] | IGS Long File Name helpers: `igs_lfn`, `sp3_filename`, `sp3_lfn`, `sp3_short_filename`, `clk_filename`, … |
//! | [`residuals_csv`] | `ResidualRecord`, `ResidualCsvWriter` (streaming) |
//! | `residuals_parquet` | `ResidualParquetWriter` (feature `parquet`) |
//! | [`manifest`] | `ManifestWriter` for [`crate::pod::run::manifest::RunManifest`] |
//! | [`qc_json`] | `QcDocument`, `write_qc_json` |
//! | [`error`] | [`PodProductsError`] |
//!
//! No estimation, propagation, or measurement modelling is performed here.
//!
//! ## References
//!
//! - Consultative Committee for Space Data Systems. (2010). Orbit Data
//!   Messages, CCSDS 502.0-B-2 / 502.0-B-3.
//! - International GNSS Service. (2020). SP3-c / SP3-d Orbit Format
//!   Specification.
#![forbid(unsafe_code)]

pub mod error;
pub mod manifest;
pub mod naming;
pub mod orbit;
pub mod qc_json;
pub mod residuals_csv;
#[cfg(feature = "pod-parquet")]
pub mod residuals_parquet;

pub use error::PodProductsError;
pub use manifest::ManifestWriter;
pub use orbit::{
    write_oem_from_spacecraft_states, write_oem_from_states, write_sp3_from_states,
    OemProductWriter, Sp3ProductWriter,
};
pub use qc_json::write_qc_json;
pub use residuals_csv::{ResidualCsvWriter, ResidualRecord};
