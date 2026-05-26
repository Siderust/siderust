// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! HTTP download with progress reporting.
//!
//! Uses `ureq` (sync, minimal dependency footprint, no async runtime).

use crate::data::{DatasetError, RuntimeDownloadMeta};
use std::io::{Read, Write};
use std::path::Path;

/// Progress callback: receives `(bytes_downloaded, total_bytes_or_zero)`.
///
/// `total` may be 0 if the server does not provide a Content-Length header.
pub type ProgressCallback = Box<dyn Fn(u64, u64)>;

/// Download a dataset to `dest`, writing through a temp file for atomicity.
///
/// `name` is used only in error messages. If `progress` is `Some`, it is
/// called periodically with `(bytes_downloaded, total_or_zero)`.
pub(super) fn download(
    name: &str,
    rdm: &RuntimeDownloadMeta,
    dest: &Path,
    progress: Option<ProgressCallback>,
) -> Result<(), DatasetError> {
    let tmp = dest.with_extension("download");

    if let Some(parent) = tmp.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let response = ureq::get(rdm.url)
        .call()
        .map_err(|e| DatasetError::Download(format!("HTTP request failed for {}: {}", name, e)))?;

    let total: u64 = response
        .headers()
        .get("Content-Length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);

    let mut reader = response.into_body().into_reader();
    let mut file = std::fs::File::create(&tmp)?;
    let mut buf = vec![0u8; 1 << 20]; // 1 MB buffer
    let mut downloaded: u64 = 0;

    loop {
        let n = reader.read(&mut buf).map_err(|e| {
            DatasetError::Download(format!("read error during {} download: {}", name, e))
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

    let file_size = std::fs::metadata(&tmp)?.len();
    if file_size < rdm.min_size {
        let _ = std::fs::remove_file(&tmp);
        return Err(DatasetError::Integrity(format!(
            "{}: downloaded file too small ({} bytes, expected >= {})",
            name, file_size, rdm.min_size,
        )));
    }

    std::fs::rename(&tmp, dest).map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        DatasetError::Io(e)
    })?;

    Ok(())
}
