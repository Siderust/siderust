// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Filesystem cache for downloaded datasets.
//!
//! Default location: `~/.siderust/data/`
//! Override with `SIDERUST_DATA_DIR` environment variable.

use super::registry::DatasetMeta;
use super::DataError;
use std::path::{Path, PathBuf};

/// Environment variable to override the data directory.
pub const DATA_DIR_ENV: &str = "SIDERUST_DATA_DIR";

/// Default subdirectory under the user's home directory.
const DEFAULT_SUBDIR: &str = ".siderust/data";

/// Resolve the data directory path.
///
/// Priority:
/// 1. `SIDERUST_DATA_DIR` environment variable
/// 2. `~/.siderust/data/`
pub fn resolve_data_dir() -> Result<PathBuf, DataError> {
    if let Ok(dir) = std::env::var(DATA_DIR_ENV) {
        let dir = dir.trim();
        if !dir.is_empty() {
            return Ok(PathBuf::from(dir));
        }
    }

    let home = std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .map_err(|_| {
            DataError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "Cannot determine home directory. Set SIDERUST_DATA_DIR explicitly.",
            ))
        })?;

    Ok(PathBuf::from(home).join(DEFAULT_SUBDIR))
}

/// Ensure the data directory exists (creates it if needed).
pub fn ensure_data_dir(dir: &Path) -> Result<(), DataError> {
    std::fs::create_dir_all(dir)?;
    Ok(())
}

/// Return the expected file path for a dataset within the cache directory.
pub fn dataset_path(data_dir: &Path, meta: &DatasetMeta) -> PathBuf {
    data_dir.join(meta.filename)
}

/// Check whether a cached dataset file exists and passes basic validation.
///
/// Returns `true` if the file exists and its size is >= `meta.min_size`.
pub fn is_cached(data_dir: &Path, meta: &DatasetMeta) -> bool {
    let path = dataset_path(data_dir, meta);
    match std::fs::metadata(&path) {
        Ok(m) => m.len() >= meta.min_size,
        Err(_) => false,
    }
}

/// Verify the integrity of a cached file.
///
/// Checks:
/// 1. File exists and size >= `min_size`
/// 2. SHA-256 matches (if `meta.sha256` is non-empty)
#[cfg(feature = "runtime-data")]
pub fn verify(path: &Path, meta: &DatasetMeta) -> Result<(), DataError> {
    let file_meta = std::fs::metadata(path)?;
    if file_meta.len() < meta.min_size {
        return Err(DataError::Integrity(format!(
            "{}: file too small ({} bytes, expected >= {})",
            meta.name,
            file_meta.len(),
            meta.min_size,
        )));
    }

    if !meta.sha256.is_empty() {
        let actual = sha256_file(path)?;
        if actual != meta.sha256 {
            return Err(DataError::Integrity(format!(
                "{}: SHA-256 mismatch (expected {}, got {})",
                meta.name, meta.sha256, actual,
            )));
        }
    }

    Ok(())
}

/// Compute SHA-256 hex digest of a file, reading in 1 MB chunks.
#[cfg(feature = "runtime-data")]
fn sha256_file(path: &Path) -> Result<String, DataError> {
    use sha2::{Digest, Sha256};
    use std::io::Read;

    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20]; // 1 MB buffer

    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }

    let digest = hasher.finalize();
    Ok(hex_encode(&digest))
}

/// Encode bytes as lowercase hex string.
#[cfg(feature = "runtime-data")]
fn hex_encode(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::registry::{DatasetId, DatasetMeta};
    use std::io::Write;

    fn dummy_meta(filename: &'static str, min_size: u64) -> DatasetMeta {
        DatasetMeta {
            id: DatasetId::De440,
            name: "test",
            url: "https://example.com/file",
            filename,
            sha256: "",
            min_size,
            size_hint: "1 B",
        }
    }

    // ── resolve_data_dir ──────────────────────────────────────────────

    #[test]
    fn resolve_data_dir_uses_env_var() {
        let tmp = std::env::temp_dir().join("siderust_test_data_dir");
        std::env::set_var(DATA_DIR_ENV, tmp.to_str().unwrap());
        let result = resolve_data_dir().unwrap();
        assert_eq!(result, tmp);
        std::env::remove_var(DATA_DIR_ENV);
    }

    #[test]
    fn resolve_data_dir_env_empty_falls_back_to_home() {
        std::env::set_var(DATA_DIR_ENV, "   ");
        let result = resolve_data_dir().unwrap();
        // Should fall back to HOME-based path since env var is whitespace-only
        let path_str = result.to_string_lossy();
        assert!(
            path_str.contains(".siderust"),
            "Expected .siderust in path: {path_str}"
        );
        std::env::remove_var(DATA_DIR_ENV);
    }

    #[test]
    fn resolve_data_dir_default_contains_siderust() {
        std::env::remove_var(DATA_DIR_ENV);
        let result = resolve_data_dir();
        // Should succeed on any machine with HOME set
        if let Ok(path) = result {
            assert!(path.to_string_lossy().contains(".siderust"));
        }
    }

    // ── ensure_data_dir ───────────────────────────────────────────────

    #[test]
    fn ensure_data_dir_creates_directory() {
        let tmp = std::env::temp_dir().join("siderust_test_ensure_dir_12345");
        let _ = std::fs::remove_dir_all(&tmp);
        ensure_data_dir(&tmp).unwrap();
        assert!(tmp.exists());
        std::fs::remove_dir_all(&tmp).ok();
    }

    #[test]
    fn ensure_data_dir_idempotent() {
        let tmp = std::env::temp_dir().join("siderust_test_ensure_dir_idempotent");
        ensure_data_dir(&tmp).unwrap();
        ensure_data_dir(&tmp).unwrap(); // should not fail on second call
        std::fs::remove_dir_all(&tmp).ok();
    }

    // ── dataset_path ──────────────────────────────────────────────────

    #[test]
    fn dataset_path_builds_correct_path() {
        let dir = std::path::PathBuf::from("/tmp/siderust_data");
        let meta = dummy_meta("de440.bsp", 0);
        let path = dataset_path(&dir, &meta);
        assert_eq!(path, dir.join("de440.bsp"));
    }

    // ── is_cached ─────────────────────────────────────────────────────

    #[test]
    fn is_cached_returns_false_when_file_absent() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_absent");
        let meta = dummy_meta("nonexistent_file.bsp", 1);
        assert!(!is_cached(&tmp, &meta));
    }

    #[test]
    fn is_cached_returns_false_when_file_too_small() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_small");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("small_file.bsp");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(b"x").unwrap();
        let meta = dummy_meta("small_file.bsp", 1000);
        assert!(!is_cached(&tmp, &meta));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    fn is_cached_returns_true_when_file_meets_min_size() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_ok");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("ok_file.bsp");
        let data: Vec<u8> = vec![0u8; 100];
        std::fs::write(&path, &data).unwrap();
        let meta = dummy_meta("ok_file.bsp", 50);
        assert!(is_cached(&tmp, &meta));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    // ── verify ────────────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "runtime-data")]
    fn verify_fails_when_file_too_small() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_small");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("tiny.bsp");
        std::fs::write(&path, b"hi").unwrap();
        let meta = dummy_meta("tiny.bsp", 1000);
        let result = verify(&path, &meta);
        assert!(result.is_err());
        let err = format!("{:?}", result.unwrap_err());
        assert!(err.contains("Integrity") || err.contains("too small"));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn verify_passes_when_no_sha256_and_size_ok() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_ok");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("ok.bsp");
        let data = vec![0u8; 200];
        std::fs::write(&path, &data).unwrap();
        let meta = dummy_meta("ok.bsp", 100);
        verify(&path, &meta).unwrap(); // should succeed with empty sha256
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn verify_fails_on_wrong_sha256() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_sha256");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("sha_check.bsp");
        let data = vec![0u8; 200];
        std::fs::write(&path, &data).unwrap();
        let meta = DatasetMeta {
            id: DatasetId::De440,
            name: "test",
            url: "https://example.com/file",
            filename: "sha_check.bsp",
            sha256: "deadbeef1234deadbeef1234deadbeef1234deadbeef1234deadbeef1234dead",
            min_size: 100,
            size_hint: "1 B",
        };
        let result = verify(&path, &meta);
        assert!(result.is_err());
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }
}
