// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Provenance metadata for bundled or computed datasets.
//!
//! Datasets ingested into a typed container (sampled spectrum, gridded
//! table, …) should carry a [`Provenance`] record describing where the
//! data came from, what version it is, and (optionally) a checksum for
//! reproducibility audits.
//!
//! Re-exported from [`crate::spectra`] and [`crate::tables`] for backwards
//! compatibility.

/// Where a dataset's samples originally came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DataSource {
    /// Cited literature (preferred for published reference data).
    LiteratureCitation {
        /// Free-form citation key (e.g. BibTeX key like `noll2012`).
        bibkey: String,
        /// Optional DOI.
        doi: Option<String>,
    },
    /// Bundled file shipped with the consuming crate.
    BundledFile {
        /// Repository-relative path of the bundled file.
        path: String,
    },
    /// External resource fetched at build- or run-time.
    External {
        /// Canonical URL of the resource.
        url: String,
    },
    /// Computed at runtime (analytic spectra, derived products, …).
    Computed {
        /// Human-readable description of the producer.
        name: String,
    },
}

/// Full provenance record for a dataset.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Provenance {
    /// Origin of the underlying samples.
    pub source: Option<DataSource>,
    /// Version, release tag, or revision string.
    pub version: Option<String>,
    /// ISO-8601 timestamp at which the data was retrieved or generated.
    pub retrieved_at: Option<String>,
    /// Optional SHA-256 checksum of the raw input as `[u8; 32]`.
    pub checksum: Option<[u8; 32]>,
    /// Free-form notes (units assumed, post-processing, caveats, …).
    pub notes: Option<String>,
}

impl Provenance {
    /// Construct an empty provenance record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Convenience constructor for a bundled file with no other metadata.
    pub fn bundled_file(path: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::BundledFile { path: path.into() }),
            ..Self::default()
        }
    }

    /// Convenience constructor for a literature citation.
    pub fn cited(bibkey: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::LiteratureCitation {
                bibkey: bibkey.into(),
                doi: None,
            }),
            ..Self::default()
        }
    }

    /// Convenience constructor for a computed dataset.
    pub fn computed(name: impl Into<String>) -> Self {
        Self {
            source: Some(DataSource::Computed { name: name.into() }),
            ..Self::default()
        }
    }

    /// Builder: attach a version string.
    pub fn with_version(mut self, version: impl Into<String>) -> Self {
        self.version = Some(version.into());
        self
    }

    /// Builder: attach free-form notes.
    pub fn with_notes(mut self, notes: impl Into<String>) -> Self {
        self.notes = Some(notes.into());
        self
    }
}
