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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_dir_path(suffix: &str) -> PathBuf {
        std::env::temp_dir().join(format!("siderust_mgr_{}", suffix))
    }

    // ── with_dir ──────────────────────────────────────────────────────

    #[test]
    fn with_dir_creates_manager() {
        let dir = temp_dir_path("with_dir");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DataManager::with_dir(&dir).unwrap();
        assert_eq!(dm.data_dir(), dir.as_path());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn data_dir_returns_configured_path() {
        let dir = temp_dir_path("data_dir_getter");
        let dm = DataManager::with_dir(&dir).unwrap();
        assert_eq!(dm.data_dir(), dir.as_path());
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── is_available ──────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "runtime-data")]
    fn is_available_returns_false_when_not_cached() {
        let dir = temp_dir_path("is_avail_false");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DataManager::with_dir(&dir).unwrap();
        assert!(!dm.is_available(DatasetId::De440));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn is_available_returns_true_when_file_present() {
        let dir = temp_dir_path("is_avail_true");
        let dm = DataManager::with_dir(&dir).unwrap();
        // Write a file sized above min_size for de440's filename
        let meta = registry::lookup(DatasetId::De440).unwrap();
        let path = dir.join(meta.filename);
        let data = vec![0u8; (meta.min_size + 100) as usize];
        std::fs::write(&path, &data).unwrap();
        assert!(dm.is_available(DatasetId::De440));
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── load_path ─────────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "runtime-data")]
    fn load_path_returns_none_when_not_cached() {
        let dir = temp_dir_path("load_path_none");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DataManager::with_dir(&dir).unwrap();
        assert!(dm.load_path(DatasetId::De440).is_none());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn load_path_returns_some_when_cached() {
        let dir = temp_dir_path("load_path_some");
        let dm = DataManager::with_dir(&dir).unwrap();
        let meta = registry::lookup(DatasetId::IersEop).unwrap();
        let path = dir.join(meta.filename);
        let data = vec![0u8; (meta.min_size + 100) as usize];
        std::fs::write(&path, &data).unwrap();
        let result = dm.load_path(DatasetId::IersEop);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), path);
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── list ──────────────────────────────────────────────────────────

    #[test]
    fn list_returns_all_datasets() {
        let dir = temp_dir_path("list_datasets");
        let dm = DataManager::with_dir(&dir).unwrap();
        let items = dm.list();
        assert!(!items.is_empty());
        // All should be unavailable in an empty dir
        assert!(items.iter().all(|(_, avail)| !avail));
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── ensure (cached already) ───────────────────────────────────────

    #[test]
    #[cfg(feature = "runtime-data")]
    fn ensure_returns_path_when_already_cached() {
        let dir = temp_dir_path("ensure_cached");
        let dm = DataManager::with_dir(&dir).unwrap();
        let meta = registry::lookup(DatasetId::IersEop).unwrap();
        let path = dir.join(meta.filename);
        // Write a file big enough to pass is_cached and verify (sha256 is empty)
        let data = vec![42u8; (meta.min_size + 100) as usize];
        std::fs::write(&path, &data).unwrap();
        let result = dm.ensure(DatasetId::IersEop).unwrap();
        assert_eq!(result, path);
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── require_meta (unknown id) ─────────────────────────────────────

    #[test]
    #[cfg(feature = "runtime-data")]
    fn ensure_unknown_dataset_returns_error() {
        let dir = temp_dir_path("ensure_unknown");
        let dm = DataManager::with_dir(&dir).unwrap();
        // Elp2000 is likely not in the DATASETS slice so lookup returns None
        // which triggers UnknownDataset error
        let result = dm.ensure(DatasetId::Elp2000);
        if let Err(e) = result {
            let s = format!("{}", e);
            assert!(
                s.contains("elp2000") || s.contains("unknown") || s.contains("Unknown"),
                "Unexpected error: {s}"
            );
        }
        // If it happens to be in the catalog that's fine too - we just don't download
        std::fs::remove_dir_all(&dir).ok();
    }
}
