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

/// Atomically write a file: write to a temporary path, then rename.
pub fn atomic_write(dest: &Path, data: &[u8]) -> Result<(), DataError> {
    let tmp = dest.with_extension("tmp");
    std::fs::write(&tmp, data)?;
    std::fs::rename(&tmp, dest)?;
    Ok(())
}
