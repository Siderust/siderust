// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Download `de440.bsp` from NAIF if not already cached.
//!
//! Uses `curl` (or `wget` as fallback) rather than adding a TLS dependency
//! to the build script. Both are available on all CI and dev machines.

use std::path::{Path, PathBuf};
use std::process::Command;

const DE440_URL: &str =
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp";
const BSP_FILENAME: &str = "de440.bsp";

/// Minimum plausible size for a valid DE440 BSP file (~114 MB).
const MIN_BSP_SIZE: u64 = 100_000_000;

/// Ensure `de440.bsp` exists in `data_dir`, downloading if necessary.
///
/// Returns the path to the BSP file.
pub fn ensure_bsp(data_dir: &Path) -> anyhow::Result<PathBuf> {
    std::fs::create_dir_all(data_dir)?;
    let bsp_path = data_dir.join(BSP_FILENAME);

    if bsp_path.exists() {
        let meta = std::fs::metadata(&bsp_path)?;
        if meta.len() > MIN_BSP_SIZE {
            eprintln!(
                "  DE440 BSP already cached ({:.1} MB)",
                meta.len() as f64 / 1_000_000.0
            );
            return Ok(bsp_path);
        }
        // Too small — probably a partial download. Re-download.
        eprintln!(
            "  DE440 BSP exists but too small ({} B), re-downloading...",
            meta.len()
        );
        std::fs::remove_file(&bsp_path)?;
    }

    eprintln!("  Downloading DE440 BSP from NAIF (~120 MB)...");
    eprintln!("  URL: {}", DE440_URL);

    // Try curl first, then wget as fallback
    let result = download_with_curl(&bsp_path)
        .or_else(|_| download_with_wget(&bsp_path));

    match result {
        Ok(()) => {
            let size = std::fs::metadata(&bsp_path)?.len();
            if size < MIN_BSP_SIZE {
                anyhow::bail!(
                    "Downloaded file too small ({} bytes), expected >{}",
                    size,
                    MIN_BSP_SIZE
                );
            }
            eprintln!("  Downloaded {:.1} MB", size as f64 / 1_000_000.0);
            Ok(bsp_path)
        }
        Err(e) => {
            anyhow::bail!(
                "Failed to download DE440 BSP. Ensure `curl` or `wget` is installed.\n\
                 You can also manually download:\n  {}\n\
                 and place it at:\n  {}\n\
                 Error: {}",
                DE440_URL,
                bsp_path.display(),
                e
            );
        }
    }
}

fn download_with_curl(dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying curl...");
    let status = Command::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "--location",           // follow redirects
            "--max-time", "900",    // 15 min timeout
            "--output",
        ])
        .arg(dest.as_os_str())
        .arg(DE440_URL)
        .status()?;

    if !status.success() {
        anyhow::bail!("curl exited with status {}", status);
    }
    Ok(())
}

fn download_with_wget(dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying wget...");
    let status = Command::new("wget")
        .args([
            "--quiet",
            "--timeout=900",
            "--output-document",
        ])
        .arg(dest.as_os_str())
        .arg(DE440_URL)
        .status()?;

    if !status.success() {
        anyhow::bail!("wget exited with status {}", status);
    }
    Ok(())
}
