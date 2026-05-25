// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level dataset manager for runtime acquisition.

use super::cache;
use super::download;
use crate::data::{lookup, Acquisition, DatasetError, DatasetId, RuntimeDownloadMeta};
use std::path::{Path, PathBuf};

/// Manages a local data cache for `RuntimeDownload` datasets.
///
/// Calling [`ensure`] or [`download`] on an `Embedded` or `ExternalProvider`
/// dataset returns [`DatasetError::NotDownloadable`].
///
/// # Example
///
/// ```rust,ignore
/// use siderust::data::{DatasetId, runtime::DatasetManager};
///
/// let dm = DatasetManager::new()?;
/// let bsp_path = dm.ensure(DatasetId::De441)?;
/// ```
pub struct DatasetManager {
    data_dir: PathBuf,
}

impl DatasetManager {
    /// Create a new `DatasetManager` using the default data directory.
    ///
    /// Default: `~/.siderust/data/` (override with `SIDERUST_DATA_DIR`).
    pub fn new() -> Result<Self, DatasetError> {
        let data_dir = cache::resolve_data_dir()?;
        cache::ensure_data_dir(&data_dir)?;
        Ok(Self { data_dir })
    }

    /// Create a `DatasetManager` pointing to a specific directory.
    pub fn with_dir(dir: impl Into<PathBuf>) -> Result<Self, DatasetError> {
        let data_dir = dir.into();
        cache::ensure_data_dir(&data_dir)?;
        Ok(Self { data_dir })
    }

    /// Returns the path to the data directory.
    pub fn data_dir(&self) -> &Path {
        &self.data_dir
    }

    /// Check whether a `RuntimeDownload` dataset is cached.
    ///
    /// Returns `false` for `Embedded` and `ExternalProvider` datasets.
    pub fn is_available(&self, id: DatasetId) -> bool {
        let meta = lookup(id);
        match &meta.acquisition {
            Acquisition::RuntimeDownload(rdm) => cache::is_cached(&self.data_dir, rdm),
            _ => false,
        }
    }

    /// Return the cached file path if the dataset is a cached `RuntimeDownload`.
    ///
    /// Returns `None` for `Embedded`, `ExternalProvider`, or uncached datasets.
    pub fn load_path(&self, id: DatasetId) -> Option<PathBuf> {
        let meta = lookup(id);
        if let Acquisition::RuntimeDownload(rdm) = &meta.acquisition {
            if cache::is_cached(&self.data_dir, rdm) {
                return Some(cache::dataset_path(&self.data_dir, rdm));
            }
        }
        None
    }

    /// Ensure a `RuntimeDownload` dataset is available: download if missing,
    /// verify integrity, and return path.
    ///
    /// Returns [`DatasetError::NotDownloadable`] for `Embedded` and
    /// `ExternalProvider` datasets.
    pub fn ensure(&self, id: DatasetId) -> Result<PathBuf, DatasetError> {
        let rdm = self.require_rdm(id)?;
        let path = cache::dataset_path(&self.data_dir, rdm);

        if cache::is_cached(&self.data_dir, rdm) {
            cache::verify(lookup(id).name, &path, rdm)?;
            return Ok(path);
        }

        self.download_inner(id, rdm, &path, None)?;
        Ok(path)
    }

    /// Download a `RuntimeDownload` dataset, with an optional progress callback.
    ///
    /// Downloads even if the file already exists (forced re-download).
    /// Returns [`DatasetError::NotDownloadable`] for other acquisition kinds.
    pub fn download(
        &self,
        id: DatasetId,
        progress: Option<download::ProgressCallback>,
    ) -> Result<PathBuf, DatasetError> {
        let rdm = self.require_rdm(id)?;
        let path = cache::dataset_path(&self.data_dir, rdm);
        self.download_inner(id, rdm, &path, progress)?;
        Ok(path)
    }

    /// List all `RuntimeDownload` datasets and their cache availability.
    pub fn list(&self) -> Vec<(DatasetId, bool)> {
        crate::data::DATASETS
            .iter()
            .filter_map(|meta| {
                if let Acquisition::RuntimeDownload(rdm) = &meta.acquisition {
                    Some((meta.id, cache::is_cached(&self.data_dir, rdm)))
                } else {
                    None
                }
            })
            .collect()
    }

    // ── Internal helpers ──────────────────────────────────────────────

    /// Extract `RuntimeDownloadMeta` or return `NotDownloadable`.
    fn require_rdm(&self, id: DatasetId) -> Result<&'static RuntimeDownloadMeta, DatasetError> {
        match &lookup(id).acquisition {
            Acquisition::RuntimeDownload(rdm) => Ok(rdm),
            _ => Err(DatasetError::NotDownloadable(id)),
        }
    }

    fn download_inner(
        &self,
        id: DatasetId,
        rdm: &'static RuntimeDownloadMeta,
        path: &Path,
        progress: Option<download::ProgressCallback>,
    ) -> Result<(), DatasetError> {
        let name = lookup(id).name;
        download::download(name, rdm, path, progress)?;
        cache::verify(name, path, rdm)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_dir_path(suffix: &str) -> PathBuf {
        std::env::temp_dir().join(format!("siderust_mgr_{}", suffix))
    }

    // ── with_dir ──────────────────────────────────────────────────────

    #[test]
    fn with_dir_creates_manager() {
        let dir = temp_dir_path("with_dir_rt");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert_eq!(dm.data_dir(), dir.as_path());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn data_dir_returns_configured_path() {
        let dir = temp_dir_path("data_dir_getter_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert_eq!(dm.data_dir(), dir.as_path());
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── is_available ──────────────────────────────────────────────────

    #[test]
    fn is_available_returns_false_when_not_cached() {
        let dir = temp_dir_path("is_avail_false_rt");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert!(!dm.is_available(DatasetId::De440));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn is_available_returns_false_for_embedded() {
        let dir = temp_dir_path("is_avail_embedded_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert!(!dm.is_available(DatasetId::Vsop87a));
        assert!(!dm.is_available(DatasetId::Elp2000));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn is_available_returns_false_for_external_provider() {
        let dir = temp_dir_path("is_avail_external_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert!(!dm.is_available(DatasetId::IersEop));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn is_available_returns_true_when_file_present() {
        let dir = temp_dir_path("is_avail_true_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let meta = lookup(DatasetId::De440);
        if let Acquisition::RuntimeDownload(rdm) = &meta.acquisition {
            let path = dir.join(rdm.filename);
            let data = vec![0u8; (rdm.min_size + 100) as usize];
            std::fs::write(&path, &data).unwrap();
            assert!(dm.is_available(DatasetId::De440));
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── load_path ─────────────────────────────────────────────────────

    #[test]
    fn load_path_returns_none_when_not_cached() {
        let dir = temp_dir_path("load_path_none_rt");
        let _ = std::fs::remove_dir_all(&dir);
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert!(dm.load_path(DatasetId::De440).is_none());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn load_path_returns_none_for_embedded() {
        let dir = temp_dir_path("load_path_embedded_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        assert!(dm.load_path(DatasetId::Vsop87a).is_none());
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── list ──────────────────────────────────────────────────────────

    #[test]
    fn list_returns_only_runtime_download_datasets() {
        let dir = temp_dir_path("list_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let items = dm.list();
        // Only De440 and De441 are RuntimeDownload
        assert_eq!(items.len(), 2);
        assert!(items.iter().any(|(id, _)| *id == DatasetId::De440));
        assert!(items.iter().any(|(id, _)| *id == DatasetId::De441));
        // All unavailable in an empty dir
        assert!(items.iter().all(|(_, avail)| !avail));
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── ensure / download on non-runtime datasets ─────────────────────

    #[test]
    fn ensure_embedded_returns_not_downloadable() {
        let dir = temp_dir_path("ensure_embedded_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let result = dm.ensure(DatasetId::Vsop87a);
        assert!(matches!(
            result,
            Err(DatasetError::NotDownloadable(DatasetId::Vsop87a))
        ));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn ensure_external_provider_returns_not_downloadable() {
        let dir = temp_dir_path("ensure_external_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let result = dm.ensure(DatasetId::IersEop);
        assert!(matches!(
            result,
            Err(DatasetError::NotDownloadable(DatasetId::IersEop))
        ));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn download_embedded_returns_not_downloadable() {
        let dir = temp_dir_path("download_embedded_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let result = dm.download(DatasetId::Elp2000, None);
        assert!(matches!(
            result,
            Err(DatasetError::NotDownloadable(DatasetId::Elp2000))
        ));
        std::fs::remove_dir_all(&dir).ok();
    }

    // ── ensure (cached already with wrong hash) ───────────────────────

    #[test]
    fn ensure_rejects_cached_file_with_wrong_pinned_hash() {
        let dir = temp_dir_path("ensure_cached_rt");
        let dm = DatasetManager::with_dir(&dir).unwrap();
        let meta = lookup(DatasetId::De440);
        if let Acquisition::RuntimeDownload(rdm) = &meta.acquisition {
            let path = dir.join(rdm.filename);
            let data = vec![42u8; (rdm.min_size + 100) as usize];
            std::fs::write(&path, &data).unwrap();
            let result = dm.ensure(DatasetId::De440);
            assert!(result.is_err());
        }
        std::fs::remove_dir_all(&dir).ok();
    }
}
