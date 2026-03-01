// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level data manager for runtime dataset acquisition.

use super::cache;
use super::download;
use super::registry::{self, DatasetId, DatasetMeta};
use super::DataError;
use std::path::{Path, PathBuf};

/// Manages a local data cache for astronomical datasets.
///
/// # Example
///
/// ```rust,ignore
/// use siderust::data::{DataManager, DatasetId};
///
/// let dm = DataManager::new()?;
/// let bsp_path = dm.ensure(DatasetId::De441)?;
/// ```
pub struct DataManager {
    data_dir: PathBuf,
}

impl DataManager {
    /// Create a new `DataManager` using the default data directory.
    ///
    /// Default: `~/.siderust/data/` (override with `SIDERUST_DATA_DIR`).
    pub fn new() -> Result<Self, DataError> {
        let data_dir = cache::resolve_data_dir()?;
        cache::ensure_data_dir(&data_dir)?;
        Ok(Self { data_dir })
    }

    /// Create a `DataManager` pointing to a specific directory.
    pub fn with_dir(dir: impl Into<PathBuf>) -> Result<Self, DataError> {
        let data_dir = dir.into();
        cache::ensure_data_dir(&data_dir)?;
        Ok(Self { data_dir })
    }

    /// Returns the path to the data directory.
    pub fn data_dir(&self) -> &Path {
        &self.data_dir
    }

    /// Check whether a dataset is already cached (no download).
    pub fn is_available(&self, id: DatasetId) -> bool {
        let Some(meta) = registry::lookup(id) else {
            return false;
        };
        cache::is_cached(&self.data_dir, meta)
    }

    /// Return the cached file path if the dataset is available, without downloading.
    pub fn load_path(&self, id: DatasetId) -> Option<PathBuf> {
        let meta = registry::lookup(id)?;
        if cache::is_cached(&self.data_dir, meta) {
            Some(cache::dataset_path(&self.data_dir, meta))
        } else {
            None
        }
    }

    /// Ensure a dataset is available: download if missing, verify integrity, return path.
    ///
    /// **This is the primary API for runtime data loading.** It will download
    /// the data only if it is not already cached — this is the user's explicit
    /// consent to download.
    pub fn ensure(&self, id: DatasetId) -> Result<PathBuf, DataError> {
        let meta = self.require_meta(id)?;
        let path = cache::dataset_path(&self.data_dir, meta);

        if cache::is_cached(&self.data_dir, meta) {
            // Already cached — verify and return
            cache::verify(&path, meta)?;
            return Ok(path);
        }

        // Download
        self.download_inner(meta, &path, None)?;
        Ok(path)
    }

    /// Download a dataset with a progress callback.
    ///
    /// Downloads even if the file already exists (forced re-download).
    pub fn download(
        &self,
        id: DatasetId,
        progress: Option<download::ProgressCallback>,
    ) -> Result<PathBuf, DataError> {
        let meta = self.require_meta(id)?;
        let path = cache::dataset_path(&self.data_dir, meta);
        self.download_inner(meta, &path, progress)?;
        Ok(path)
    }

    /// List all datasets and their cache status.
    pub fn list(&self) -> Vec<(DatasetId, bool)> {
        registry::DATASETS
            .iter()
            .map(|meta| (meta.id, cache::is_cached(&self.data_dir, meta)))
            .collect()
    }

    // ── Internal helpers ──────────────────────────────────────────────

    fn require_meta(&self, id: DatasetId) -> Result<&'static DatasetMeta, DataError> {
        registry::lookup(id).ok_or_else(|| DataError::UnknownDataset(id.as_str().to_string()))
    }

    fn download_inner(
        &self,
        meta: &DatasetMeta,
        path: &Path,
        progress: Option<download::ProgressCallback>,
    ) -> Result<(), DataError> {
        download::download(meta, path, progress)?;
        cache::verify(path, meta)?;
        Ok(())
    }
}
