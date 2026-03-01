// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! HTTP download with progress reporting.
//!
//! Uses `ureq` (sync, minimal dependency footprint — no async runtime).

use super::registry::DatasetMeta;
use super::DataError;
use std::io::{Read, Write};
use std::path::Path;

/// Progress callback: receives `(bytes_downloaded, total_bytes_or_zero)`.
///
/// `total` may be 0 if the server does not provide a Content-Length header.
pub type ProgressCallback = Box<dyn Fn(u64, u64)>;

/// Download a dataset to `dest`, writing through a temp file for atomicity.
///
/// If `progress` is `Some`, it is called periodically with download progress.
pub fn download(
    meta: &DatasetMeta,
    dest: &Path,
    progress: Option<ProgressCallback>,
) -> Result<(), DataError> {
    let tmp = dest.with_extension("download");

    // Ensure parent directory exists
    if let Some(parent) = tmp.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let response = ureq::get(meta.url).call().map_err(|e| {
        DataError::Download(format!("HTTP request failed for {}: {}", meta.name, e))
    })?;

    let total: u64 = response
        .header("Content-Length")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);

    let mut reader = response.into_reader();
    let mut file = std::fs::File::create(&tmp)?;
    let mut buf = vec![0u8; 1 << 20]; // 1 MB buffer
    let mut downloaded: u64 = 0;

    loop {
        let n = reader.read(&mut buf).map_err(|e| {
            DataError::Download(format!("read error during {} download: {}", meta.name, e))
        })?;
        if n == 0 {
            break;
        }
        file.write_all(&buf[..n])?;
        downloaded += n as u64;

        if let Some(ref cb) = progress {
            cb(downloaded, total);
        }
    }

    file.flush()?;
    drop(file);

    // Basic size check before renaming into place
    let file_size = std::fs::metadata(&tmp)?.len();
    if file_size < meta.min_size {
        let _ = std::fs::remove_file(&tmp);
        return Err(DataError::Integrity(format!(
            "{}: downloaded file too small ({} bytes, expected >= {})",
            meta.name, file_size, meta.min_size,
        )));
    }

    // Atomic rename
    std::fs::rename(&tmp, dest).map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        DataError::Io(e)
    })?;

    Ok(())
}
