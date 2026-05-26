// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Provenance metadata
//!
//! ## Scientific scope
//!
//! Reproducible astronomical pipelines must be able to answer the
//! question "where did this number come from?" for every dataset they
//! consume — was it the literature curve cited in a paper, a specific
//! release of a vendored ASCII table, the output of an external service
//! pulled at build time, or a derived product computed at runtime?
//! The provenance record carried by every typed dataset in this crate
//! lets downstream consumers audit the literature source, file version,
//! and (optionally) cryptographic hash of the bytes that produced any
//! computed result, without having to re-derive that information from
//! file paths or imported strings.
//!
//! This is the same hygiene principle followed by the IVOA Provenance
//! Data Model and by IERS conventions when distributing reference data.
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - [`Provenance`] — full record (origin, version, retrieval timestamp,
//!   optional SHA-256 checksum, free-form notes).
//! - [`DataSource`] — origin classifier with variants for
//!   `LiteratureCitation`, `BundledFile`, `External`, and `Computed`.
//! - Builder helpers: [`Provenance::new`], [`Provenance::bundled_file`],
//!   [`Provenance::cited`], [`Provenance::computed`],
//!   [`Provenance::with_version`], [`Provenance::with_notes`].
//!
//! All fields are owned `String` / `Option<String>` for portability; no
//! external serialization is implied. Re-exported at the crate root and
//! from the `photometry` and `tables` feature modules for backwards
//! compatibility.
//!
//! The [`checksum`] submodule provides a const-evaluable SHA-256 and the
//! [`assert_data_checksum!`](crate::assert_data_checksum) macro for
//! pinning the hash of any [`include_str!`] / [`include_bytes!`] data
//! blob shipped inside the crate. Mismatches become hard compile errors.
//!
//! ## References
//!
//! - International Virtual Observatory Alliance (2020). *IVOA Provenance
//!   Data Model*. IVOA Recommendation, version 1.0.
//!   <https://www.ivoa.net/documents/ProvenanceDM/>.
//! - IERS Conventions (2010). *IERS Technical Note 36*, Chapter 1
//!   (citing reference data products).

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
