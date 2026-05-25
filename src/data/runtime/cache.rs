// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Filesystem cache for downloaded datasets.
//!
//! Default location: `~/.siderust/data/`
//! Override with `SIDERUST_DATA_DIR` environment variable.

use crate::data::DatasetError;
use crate::data::RuntimeDownloadMeta;
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
pub fn resolve_data_dir() -> Result<PathBuf, DatasetError> {
    if let Ok(dir) = std::env::var(DATA_DIR_ENV) {
        let dir = dir.trim();
        if !dir.is_empty() {
            return Ok(PathBuf::from(dir));
        }
    }

    let home = std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .map_err(|_| {
            DatasetError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "Cannot determine home directory. Set SIDERUST_DATA_DIR explicitly.",
            ))
        })?;

    Ok(PathBuf::from(home).join(DEFAULT_SUBDIR))
}

/// Ensure the data directory exists (creates it if needed).
pub fn ensure_data_dir(dir: &Path) -> Result<(), DatasetError> {
    std::fs::create_dir_all(dir)?;
    Ok(())
}

/// Return the expected file path for a dataset within the cache directory.
pub fn dataset_path(data_dir: &Path, rdm: &RuntimeDownloadMeta) -> PathBuf {
    data_dir.join(rdm.filename)
}

/// Check whether a cached dataset file exists and passes basic size validation.
///
/// Returns `true` if the file exists and its size is >= [`rdm.min_size`].
pub fn is_cached(data_dir: &Path, rdm: &RuntimeDownloadMeta) -> bool {
    let path = dataset_path(data_dir, rdm);
    match std::fs::metadata(&path) {
        Ok(m) => m.len() >= rdm.min_size,
        Err(_) => false,
    }
}

/// Verify the integrity of a cached file.
///
/// Checks:
/// 1. File exists and size >= `min_size`
/// 2. SHA-256 matches (if `rdm.sha256` is non-empty)
pub fn verify(name: &str, path: &Path, rdm: &RuntimeDownloadMeta) -> Result<(), DatasetError> {
    let file_meta = std::fs::metadata(path)?;
    if file_meta.len() < rdm.min_size {
        return Err(DatasetError::Integrity(format!(
            "{}: file too small ({} bytes, expected >= {})",
            name,
            file_meta.len(),
            rdm.min_size,
        )));
    }

    if !rdm.sha256.is_empty() {
        let actual = sha256_file(path)?;
        if actual != rdm.sha256 {
            return Err(DatasetError::Integrity(format!(
                "{}: SHA-256 mismatch (expected {}, got {})",
                name, rdm.sha256, actual,
            )));
        }
    }

    Ok(())
}

/// Compute SHA-256 hex digest of a file, reading in 1 MB chunks.
fn sha256_file(path: &Path) -> Result<String, DatasetError> {
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
fn hex_encode(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn dummy_rdm(filename: &'static str, min_size: u64) -> RuntimeDownloadMeta {
        RuntimeDownloadMeta {
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
        if let Ok(path) = result {
            assert!(path.to_string_lossy().contains(".siderust"));
        }
    }

    // ── ensure_data_dir ───────────────────────────────────────────────

    #[test]
    fn ensure_data_dir_creates_directory() {
        let tmp = std::env::temp_dir().join("siderust_test_ensure_dir_12345_rt");
        let _ = std::fs::remove_dir_all(&tmp);
        ensure_data_dir(&tmp).unwrap();
        assert!(tmp.exists());
        std::fs::remove_dir_all(&tmp).ok();
    }

    #[test]
    fn ensure_data_dir_idempotent() {
        let tmp = std::env::temp_dir().join("siderust_test_ensure_dir_idempotent_rt");
        ensure_data_dir(&tmp).unwrap();
        ensure_data_dir(&tmp).unwrap();
        std::fs::remove_dir_all(&tmp).ok();
    }

    // ── dataset_path ──────────────────────────────────────────────────

    #[test]
    fn dataset_path_builds_correct_path() {
        let dir = std::path::PathBuf::from("/tmp/siderust_data");
        let rdm = dummy_rdm("de440.bsp", 0);
        let path = dataset_path(&dir, &rdm);
        assert_eq!(path, dir.join("de440.bsp"));
    }

    // ── is_cached ─────────────────────────────────────────────────────

    #[test]
    fn is_cached_returns_false_when_file_absent() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_absent_rt");
        let rdm = dummy_rdm("nonexistent_file.bsp", 1);
        assert!(!is_cached(&tmp, &rdm));
    }

    #[test]
    fn is_cached_returns_false_when_file_too_small() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_small_rt");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("small_file.bsp");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(b"x").unwrap();
        let rdm = dummy_rdm("small_file.bsp", 1000);
        assert!(!is_cached(&tmp, &rdm));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    fn is_cached_returns_true_when_file_meets_min_size() {
        let tmp = std::env::temp_dir().join("siderust_test_cache_ok_rt");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("ok_file.bsp");
        let data: Vec<u8> = vec![0u8; 100];
        std::fs::write(&path, &data).unwrap();
        let rdm = dummy_rdm("ok_file.bsp", 50);
        assert!(is_cached(&tmp, &rdm));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    // ── verify ────────────────────────────────────────────────────────

    #[test]
    fn verify_fails_when_file_too_small() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_small_rt");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("tiny.bsp");
        std::fs::write(&path, b"hi").unwrap();
        let rdm = dummy_rdm("tiny.bsp", 1000);
        let result = verify("test", &path, &rdm);
        assert!(result.is_err());
        let err = format!("{:?}", result.unwrap_err());
        assert!(err.contains("Integrity") || err.contains("too small"));
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    fn verify_passes_when_no_sha256_and_size_ok() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_ok_rt");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("ok.bsp");
        let data = vec![0u8; 200];
        std::fs::write(&path, &data).unwrap();
        let rdm = dummy_rdm("ok.bsp", 100);
        verify("test", &path, &rdm).unwrap();
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }

    #[test]
    fn verify_fails_on_wrong_sha256() {
        let tmp = std::env::temp_dir().join("siderust_test_verify_sha256_rt");
        ensure_data_dir(&tmp).unwrap();
        let path = tmp.join("sha_check.bsp");
        let data = vec![0u8; 200];
        std::fs::write(&path, &data).unwrap();
        let rdm = RuntimeDownloadMeta {
            url: "https://example.com/file",
            filename: "sha_check.bsp",
            sha256: "deadbeef1234deadbeef1234deadbeef1234deadbeef1234deadbeef1234dead",
            min_size: 100,
            size_hint: "1 B",
        };
        let result = verify("test", &path, &rdm);
        assert!(result.is_err());
        std::fs::remove_file(&path).ok();
        std::fs::remove_dir(&tmp).ok();
    }
}
