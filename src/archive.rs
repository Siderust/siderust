// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Archive registry access layer
//!
//! This module exposes the build-time generated [`ARCHIVE_ENTRIES`] table that
//! mirrors the contents of `archive/MANIFEST.toml`. The registry is populated
//! by `build.rs` only when the `archive-data` feature is enabled **and** the
//! `archive/` submodule is checked out. Otherwise it is an empty slice.
//!
//! Consumers use this module to enumerate available datasets, look them up by
//! family id, and resolve the on-disk paths that point into the archive.
//!
//! ## Feature gating
//!
//! * `archive-data` — enables manifest parsing and populates the registry.
//! * `embedded-data` — extends `archive-data` with `include_bytes!` payloads
//!   for selected datasets.
//! * `external-data` — exposes [`AchiveLocator`] helpers for loading datasets
//!   from a user-specified base directory.
//!
//! Without any of those features the module still exists, [`ARCHIVE_ENTRIES`]
//! is `&[]`, and [`lookup_family`] returns `None`.

use core::fmt;

/// A single archive registry entry, generated from the per-family manifests.
#[derive(Debug, Clone, Copy)]
pub struct ArchiveEntry {
    /// Family identifier (matches `[[family]].id` in `archive/MANIFEST.toml`).
    pub family: &'static str,
    /// Dataset identifier (matches `dataset_id` in the family manifest).
    pub dataset_id: &'static str,
    /// Dataset kind (e.g. `lagrange-chebyshev`, `planetary-theory`).
    pub kind: &'static str,
    /// Path of the family directory relative to the archive root.
    pub relative_path: &'static str,
    /// Path of the family manifest relative to the archive root.
    pub manifest_path: &'static str,
    /// Inclusive lower bound of the validity interval in Julian Date (UT-agnostic).
    pub valid_from_jd: f64,
    /// Inclusive upper bound of the validity interval in Julian Date.
    pub valid_to_jd: f64,
    /// Aggregate SHA-256 checksum of the family payload (may be empty).
    pub checksum_sha256: &'static str,
}

impl fmt::Display for ArchiveEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{family} ({kind}) → {path}",
            family = self.family,
            kind = self.kind,
            path = self.relative_path,
        )
    }
}

include!(concat!(env!("OUT_DIR"), "/archive_registry.rs"));

/// Returns the registry entry for the supplied family id, if present.
#[must_use]
pub fn lookup_family(family: &str) -> Option<&'static ArchiveEntry> {
    ARCHIVE_ENTRIES.iter().find(|entry| entry.family == family)
}

/// Returns `true` if the archive registry was populated at build time.
#[must_use]
pub fn registry_populated() -> bool {
    !ARCHIVE_ENTRIES.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_consistent() {
        // The registry is always present (possibly empty) and lookups
        // never panic regardless of feature configuration.
        let entries = ARCHIVE_ENTRIES;
        for entry in entries {
            assert!(!entry.family.is_empty());
            assert!(!entry.kind.is_empty());
        }
        assert!(lookup_family("__missing__").is_none());
    }
}
