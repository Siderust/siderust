// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! VSOP87 **dataset fetcher**
//!
//! First stage of the build pipeline: make sure the raw data files are present
//! locally.  If *any* VSOP87 file is already found under `data_dir/` the stage
//! is considered satisfied; otherwise all files are downloaded from the IMCCE
//! public FTP mirror (served over HTTPS as well).
//!
//! ```text
//! Remote  : https://ftp.imcce.fr/pub/ephem/planets/vsop87/
//! Local   : <workspace>/dataset/
//! Trigger : called by `vsop87::run()` **before** `collect_terms()`.
//! ```
//!
//! ## Crate features required in `Cargo.toml`
//! Add under **[build-dependencies]** because this code runs inside the build
//! script:
//! ```toml
//! reqwest   = { version = "0.12", features = ["blocking", "rustls-tls"] }
//! regex     = "1"
//! anyhow    = "1"
//! ```
//! The `blocking` client is easier inside build scripts (no async runtime).
//!
//! ## High‑level algorithm
//! 1. **Quick check**: if `data_dir` already contains at least one
//!    `VSOP87X.xxx` file → return early.
//! 2. **Local copy**: attempt to copy the dataset from `scripts/vsop87/dataset`
//!    (checked in via Git LFS).
//! 3. **Fetch directory listing** from the host and extract all matching file
//!    names via regex (only if local copy was missing).
//! 4. **Download** each file only if it does not yet exist (makes rebuilds
//!    faster).
//!
//! Any network or filesystem error bubbles up as `anyhow::Error` so the build
//! fails loudly instead of producing incomplete tables.

use std::{
    env,
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use regex::Regex;

const BASE_URL: &str = "https://ftp.imcce.fr/pub/ephem/planets/vsop87/";
const FILE_RE_STR: &str = r"VSOP87[A-Z]\.[A-Za-z]{3}"; // e.g. VSOP87A.ear

/// Public entry: ensure that `data_dir` contains the VSOP87 dataset.
///
/// * Creates the directory if it doesn’t exist.
/// * Tries to copy files from the repo’s `dataset/` folder (Git LFS).
/// * Downloads any missing file from the internet as a fallback.
pub fn ensure_dataset(data_dir: &Path) -> Result<()> {
    // 1) Early‑out if dataset already present -----------------------------
    if contains_vsop_files(data_dir)? {
        return Ok(());
    }

    // 2) Try local checkout (Git LFS) ------------------------------------
    if let Ok(true) = copy_from_repo(data_dir) {
        return Ok(());
    }

    fs::create_dir_all(data_dir)
        .with_context(|| format!("Could not create directory {data_dir:?}"))?;

    // 3) Obtain remote file list ------------------------------------------
    let files = remote_file_list().context("Fetching VSOP87 directory listing")?;

    // 4) Download each file (skipping those already present) --------------
    for file in files {
        let dst = data_dir.join(&file);
        if dst.exists() {
            continue; // keep local copy – useful during incremental builds
        }
        let url = format!("{BASE_URL}{file}");
        download_file(&url, &dst)?;
    }
    Ok(())
}

/// Copy dataset files from the repository if available.
///
/// Returns `Ok(true)` if files were copied successfully.
fn copy_from_repo(dst: &Path) -> Result<bool> {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/vsop87/dataset");
    eprintln!("Attempting to copy VSOP87 data from: {:?}", src);
    eprintln!("Target directory: {:?}", dst);
    eprintln!("Source exists: {}", src.exists());

    if !contains_vsop_files(&src)? {
        eprintln!("No VSOP87 files found in source directory");
        return Ok(false);
    }

    eprintln!("Found VSOP87 files in source, copying...");
    fs::create_dir_all(dst).with_context(|| format!("Could not create directory {dst:?}"))?;

    let mut copied_count = 0;
    for entry in fs::read_dir(&src)? {
        let entry = entry?;
        let target = dst.join(entry.file_name());
        fs::copy(entry.path(), target)
            .with_context(|| format!("copy {:?} -> {:?}", entry.path(), dst))?;
        copied_count += 1;
    }
    eprintln!("Successfully copied {} files", copied_count);
    Ok(true)
}

/// Check quickly whether the directory already contains at least one VSOP87
/// file.  Saving a network round‑trip in the common incremental‑build case.
fn contains_vsop_files(dir: &Path) -> Result<bool> {
    if !dir.is_dir() {
        return Ok(false);
    }
    let re = Regex::new(FILE_RE_STR).unwrap();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if re.is_match(&entry.file_name().to_string_lossy()) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Retrieve the directory listing from the IMCCE server and return every file
/// name that matches [`FILE_RE_STR`].
fn remote_file_list() -> Result<Vec<String>> {
    let body = reqwest::blocking::get(BASE_URL)?.text()?;
    let re = Regex::new(FILE_RE_STR).unwrap();
    let mut files: Vec<String> = re.find_iter(&body).map(|m| m.as_str().to_owned()).collect();
    files.sort();
    files.dedup();
    Ok(files)
}

/// Download a single file `url` into path `dst`.
fn download_file(url: &str, dst: &PathBuf) -> Result<()> {
    println!("cargo:info=Downloading {} → {}", url, dst.display());
    let bytes = reqwest::blocking::get(url)?.bytes()?;
    let mut fh = File::create(dst).with_context(|| format!("Failed to create {dst:?}"))?;
    fh.write_all(&bytes)
        .with_context(|| format!("Failed to write {dst:?}"))?;
    Ok(())
}
