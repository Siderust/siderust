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
