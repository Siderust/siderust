// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! VSOP87 **output stage**
//!
//! `io.rs` is responsible only for **side‑effects**: taking the Rust source
//! blobs produced by [`codegen`](crate::scripts::vsop87::codegen) and writing
//! one file per VSOP87 version (`vsop87a.rs`, `vsop87e.rs`, …) inside the
//! directory pointed at by Cargo’s `OUT_DIR` environment variable.
//!
//! Why a separate module?
//! * Keeps pure/functional stages (`collect`, `codegen`) free of I/O.
//! * Makes it easy to swap the backend later (e.g. in‑memory, mock FS, etc.).
//!
//! The function exposed here is minimal: it receives the map returned by
//! `codegen::generate_modules` and the already‑resolved `out_dir` path; it
//! performs the write and surfaces any error through `anyhow`.

use std::{collections::BTreeMap, fs::File, io::Write, path::Path};

use anyhow::Context;

/// Persist each `(version → source_code)` pair into an actual `.rs` file.
///
/// * `mods` – map keyed by the **version letter** (A, E, …) produced by
///   `codegen`.
/// * `out_dir` – directory where Cargo expects build scripts to drop their
///   generated artefacts.  Usually passed from the caller as
///   `PathBuf::from(env!("OUT_DIR"))`.
///
/// For every `(version, code)` entry we create a file named
/// `vsop87{version}.rs`, lowering the version letter for consistency, e.g.:
///
/// ```text
/// vsop87a.rs  // ← version 'A'
/// vsop87e.rs  // ← version 'E'
/// ```
///
/// The function is **idempotent**: it always overwrites the file with the new
/// contents.  Errors from the file system are wrapped with context so that the
/// caller knows *which* file failed.
pub fn write_modules(mods: &BTreeMap<char, String>, out_dir: &Path) -> anyhow::Result<()> {
    for (version, code) in mods {
        // ------------------------------------------------------------------
        // 1) Compute final path   OUT_DIR / "vsop87{version}.rs"
        // ------------------------------------------------------------------
        let file_name = format!("vsop87{}.rs", version.to_ascii_lowercase());
        let path = out_dir.join(&file_name);

        // ------------------------------------------------------------------
        // 2) Create (or truncate) the file and write the source string
        // ------------------------------------------------------------------
        let mut f = File::create(&path).with_context(|| format!("Could not create {path:?}"))?;
        f.write_all(code.as_bytes())
            .with_context(|| format!("Error writing {path:?}"))?;

        // A friendly note for `cargo build -vv` users.
        println!("cargo:info=Generated {file_name}");
    }
    Ok(())
}
