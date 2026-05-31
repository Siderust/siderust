// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Runtime data management
//!
//! ## Scientific scope
//!
//! Several reference datasets used in modern astronomy — JPL planetary
//! ephemerides (DE440, DE441), high-resolution IERS time series, large
//! star catalogues — are tens to hundreds of megabytes and cannot be
//! shipped inside a Rust crate without inflating every dependent build
//! and bloating the registry. The community-standard solution, used by
//! Astropy and SPICE, is to keep the heavy data *out* of the source
//! distribution and instead download it on demand to a per-user cache,
//! verifying integrity by SHA-256 against a pinned manifest.
//!
//! This module implements that pattern for `siderust`. It does not by
//! itself contain any astronomical algorithm; rather, it materialises the
//! coefficient files that ephemeris and IERS modules consume, so that
//! results are byte-identical to what one would obtain from the
//! upstream Astropy / NAIF SPICE distributions.
//!
//! ## Technical scope
//!
//! Behind the `runtime-data` feature this module provides:
//!
//! - `DataManager` — discovers (or creates) the cache directory,
//!   downloads, verifies, and returns paths to dataset files.
//! - [`DatasetId`], [`DatasetMeta`], [`DATASETS`] — pinned manifest of
//!   the datasets known to the crate.
//! - [`DataError`] — error taxonomy (`Io`, `Download`, `Integrity`,
//!   `UnknownDataset`, `Parse`).
//! - `ProgressCallback` — optional progress hook for downloads.
//!
//! Without the `runtime-data` feature, only [`DatasetId`],
//! [`DatasetMeta`], [`DATASETS`], and [`DataError`] are available, so
//! that downstream code can still reason about the *registry* without
//! pulling in the network/IO stack.
//!
//! Cache location is `~/.siderust/data/` by default; override with the
//! `SIDERUST_DATA_DIR` environment variable. Data is never downloaded
//! without an explicit call to `DataManager::ensure` /
//! `DataManager::download`.
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use siderust::data::{DataManager, DatasetId};
//!
//! let dm = DataManager::new()?;
//! let path = dm.ensure(DatasetId::De441)?;
//! println!("DE441 ready at: {}", path.display());
//! ```
//!
//! ## References
//!
//! - The Astropy Collaboration (2022). "The Astropy Project: Sustaining
//!   and Growing a Community-oriented Open-Source Project and the Latest
//!   Major Release (v5.0) of the Core Package". *The Astrophysical
//!   Journal* **935**, 167. doi:10.3847/1538-4357/ac7c74.
//! - Acton, C. H., Bachman, N., Semenov, B., Wright, E. (2018). "A look
//!   towards the future in the handling of space science mission
//!   geometry". *Planetary and Space Science* **150**, 9–12.
//!   doi:10.1016/j.pss.2017.02.013.
//! - NAIF (2022). *SPK Required Reading* and *DAF Required Reading*.
//!   <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>.

#[cfg(feature = "runtime-data")]
mod cache;
mod registry;

#[cfg(feature = "runtime-data")]
mod download;

pub mod daf;
pub mod spk;
mod spk_kernel;

#[cfg(feature = "runtime-data")]
mod manager;

pub use registry::{DatasetId, DatasetMeta, DATASETS};
pub use spk_kernel::{SpkKernelError, SpkKernelSet};

#[cfg(feature = "runtime-data")]
pub use download::ProgressCallback;
#[cfg(feature = "runtime-data")]
pub use manager::DataManager;

/// Errors that can occur during data operations.
#[derive(Debug)]
pub enum DataError {
    /// I/O error (file system, permissions, etc.).
    Io(std::io::Error),
    /// HTTP download error.
    #[cfg(feature = "runtime-data")]
    Download(String),
    /// Data integrity check failed (checksum mismatch, wrong size, etc.).
    Integrity(String),
    /// The requested dataset was not found in the registry.
    UnknownDataset(String),
    /// DAF/SPK parsing error.
    Parse(String),
}

impl std::fmt::Display for DataError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DataError::Io(e) => write!(f, "I/O error: {}", e),
            #[cfg(feature = "runtime-data")]
            DataError::Download(msg) => write!(f, "download error: {}", msg),
            DataError::Integrity(msg) => write!(f, "integrity check failed: {}", msg),
            DataError::UnknownDataset(id) => write!(f, "unknown dataset: {}", id),
            DataError::Parse(msg) => write!(f, "parse error: {}", msg),
        }
    }
}

impl std::error::Error for DataError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DataError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for DataError {
    fn from(e: std::io::Error) -> Self {
        DataError::Io(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_error_display_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let e = DataError::Io(io_err);
        let msg = format!("{}", e);
        assert!(msg.contains("I/O error"));
    }

    #[test]
    fn data_error_display_integrity() {
        let e = DataError::Integrity("checksum mismatch".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("integrity check failed"));
        assert!(msg.contains("checksum mismatch"));
    }

    #[test]
    fn data_error_display_unknown_dataset() {
        let e = DataError::UnknownDataset("de999".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("unknown dataset"));
        assert!(msg.contains("de999"));
    }

    #[test]
    fn data_error_display_parse() {
        let e = DataError::Parse("bad header".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("parse error"));
        assert!(msg.contains("bad header"));
    }

    #[test]
    fn data_error_from_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "denied");
        let e: DataError = io_err.into();
        matches!(e, DataError::Io(_));
    }

    #[test]
    fn data_error_source_io_is_some() {
        use std::error::Error;
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "not found");
        let e = DataError::Io(io_err);
        assert!(e.source().is_some());
    }

    #[test]
    fn data_error_source_non_io_is_none() {
        use std::error::Error;
        let e = DataError::Integrity("x".to_string());
        assert!(e.source().is_none());
    }

    #[test]
    fn data_error_debug() {
        let e = DataError::Parse("test".to_string());
        let dbg = format!("{:?}", e);
        assert!(dbg.contains("Parse"));
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn data_error_display_download() {
        let e = DataError::Download("timeout".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("download error"));
        assert!(msg.contains("timeout"));
    }
}
