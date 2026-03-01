// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Runtime Data Management
//!
//! This module provides runtime download, caching, and loading of astronomical
//! datasets. It follows an approach similar to [Astropy's data management][astropy]:
//! heavy datasets are **not** embedded at compile time but downloaded on demand
//! to a local cache directory (`~/.siderust/data/` by default).
//!
//! [astropy]: https://docs.astropy.org/en/stable/utils/data.html
//!
//! ## Design Principles
//!
//! - **Explicit consent**: Data is never downloaded without the user calling
//!   [`DataManager::ensure`] or [`DataManager::download`] explicitly.
//! - **Persistent cache**: Downloaded files are stored on disk and reused
//!   across program runs.
//! - **Integrity verification**: Downloads are validated via SHA-256 checksums
//!   and minimum file sizes.
//! - **Environment override**: Set `SIDERUST_DATA_DIR` to relocate the cache.
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use siderust::data::{DataManager, DatasetId};
//!
//! let dm = DataManager::new()?;
//!
//! // Download DE441 BSP if not already cached (explicit consent)
//! let path = dm.ensure(DatasetId::De441)?;
//! println!("DE441 ready at: {}", path.display());
//! ```
//!
//! ## Feature Gate
//!
//! This module requires the `runtime-data` Cargo feature:
//! ```toml
//! siderust = { version = "0.5", features = ["runtime-data"] }
//! ```

#[cfg(feature = "runtime-data")]
mod cache;
mod registry;

#[cfg(feature = "runtime-data")]
mod download;

pub mod daf;
pub mod spk;

#[cfg(feature = "runtime-data")]
mod manager;

pub use registry::{DatasetId, DatasetMeta, DATASETS};

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
