// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Download `de441_part-2.bsp` from NAIF if not already cached.
//!
//! Uses `curl` (or `wget` as fallback) rather than adding a TLS dependency
//! to the build script. Both are available on all CI and dev machines.

use std::path::{Path, PathBuf};
use std::process::Command;

const DE441_URL: &str =
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp";
const BSP_FILENAME: &str = "de441_part-2.bsp";

/// Minimum plausible size for a valid DE441 part-2 BSP file (~1.6 GB).
const MIN_BSP_SIZE: u64 = 1_500_000_000;

/// Ensure `de441_part-2.bsp` exists in `data_dir`, downloading if necessary.
///
/// Returns the path to the BSP file.
///
/// Algorithm:
/// 1. If file already exists in data_dir with valid size → return early
/// 2. Try to copy from `scripts/de441/dataset/` (Git LFS backup)
/// 3. Download from NAIF as fallback
pub fn ensure_bsp(data_dir: &Path) -> anyhow::Result<PathBuf> {
    std::fs::create_dir_all(data_dir)?;
    let bsp_path = data_dir.join(BSP_FILENAME);

    // 1) Early return if already present with valid size
    if bsp_path.exists() {
        let meta = std::fs::metadata(&bsp_path)?;
        if meta.len() > MIN_BSP_SIZE {
            eprintln!(
                "  DE441 BSP already cached ({:.1} MB)",
                meta.len() as f64 / 1_000_000.0
            );
            return Ok(bsp_path);
        }
        // Too small — probably a partial download. Remove and re-acquire.
        eprintln!(
            "  DE441 BSP exists but too small ({} B), re-acquiring...",
            meta.len()
        );
        std::fs::remove_file(&bsp_path)?;
    }

    // 2) Try local checkout (Git LFS)
    if let Ok(true) = copy_from_repo(&bsp_path) {
        return Ok(bsp_path);
    }

    // 3) Download from NAIF as fallback
    eprintln!("  Downloading DE441 part-2 BSP from NAIF (~1.65 GB)...");
    eprintln!("  URL: {}", DE441_URL);

    // Try curl first, then wget as fallback
    let result = download_with_curl(&bsp_path).or_else(|_| download_with_wget(&bsp_path));

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
                "Failed to download DE441 part-2 BSP. Ensure `curl` or `wget` is installed.\n\
                 You can also manually download:\n  {}\n\
                 and place it at:\n  {}\n\
                 Error: {}",
                DE441_URL,
                bsp_path.display(),
                e
            );
        }
    }
}

/// Copy de441_part-2.bsp from the repository if available.
///
/// Returns `Ok(true)` if the file was copied successfully with valid size.
fn copy_from_repo(dst: &Path) -> anyhow::Result<bool> {
    let src = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("scripts/de441/dataset")
        .join(BSP_FILENAME);

    eprintln!(
        "  Attempting to copy DE441 part-2 BSP from Git LFS: {:?}",
        src
    );

    if !src.exists() {
        eprintln!("  No local BSP found in source tree");
        return Ok(false);
    }

    let meta = std::fs::metadata(&src)?;
    if meta.len() < MIN_BSP_SIZE {
        eprintln!("  Local BSP file too small ({} B), skipping", meta.len());
        return Ok(false);
    }

    eprintln!(
        "  Found valid BSP in source ({:.1} MB), copying...",
        meta.len() as f64 / 1_000_000.0
    );
    std::fs::copy(&src, dst)?;
    eprintln!("  Successfully copied from Git LFS");
    Ok(true)
}

fn download_with_curl(dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying curl...");
    let status = Command::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "--location", // follow redirects
            "--max-time",
            "5400", // 90 min timeout
            "--output",
        ])
        .arg(dest.as_os_str())
        .arg(DE441_URL)
        .status()?;

    if !status.success() {
        anyhow::bail!("curl exited with status {}", status);
    }
    Ok(())
}

fn download_with_wget(dest: &Path) -> anyhow::Result<()> {
    eprintln!("  Trying wget...");
    let status = Command::new("wget")
        .args(["--quiet", "--timeout=5400", "--output-document"])
        .arg(dest.as_os_str())
        .arg(DE441_URL)
        .status()?;

    if !status.success() {
        anyhow::bail!("wget exited with status {}", status);
    }
    Ok(())
}
