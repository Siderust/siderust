// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Scientific dataset catalog and runtime acquisition
//!
//! This module is the umbrella catalog for all scientific datasets known to
//! `siderust`. It answers two questions:
//!
//! 1. **What datasets exist?** — [`DatasetId`], [`DatasetMeta`], [`DATASETS`],
//!    [`lookup`].
//! 2. **How is each dataset acquired?** — [`Acquisition`] discriminates between
//!    datasets that are compiled-in ([`Acquisition::Embedded`]), those that must
//!    be downloaded at runtime
//!    ([`Acquisition::RuntimeDownload`]), and those owned by another crate
//!    ([`Acquisition::ExternalProvider`]).
//!
//! ## Runtime data management
//!
//! Behind the `runtime-data` feature the [`runtime`] sub-module exposes
//! [`runtime::DatasetManager`] for downloading and caching
//! `RuntimeDownload` datasets. Calling acquire methods on `Embedded` or
//! `ExternalProvider` datasets returns [`DatasetError::NotDownloadable`].
//!
//! ## Error handling
//!
//! All dataset-level errors are reported as [`DatasetError`]. Parse errors from
//! the SPICE format layer arrive wrapped in [`DatasetError::Spice`] via the
//! [`From<formats::spice::SpiceError>`] conversion.
//!
//! ## Quick start
//!
//! ```rust,ignore
//! use siderust::datasets::{DatasetId, runtime::DatasetManager};
//!
//! let dm = DatasetManager::new()?;
//! let path = dm.ensure(DatasetId::De441)?;
//! println!("DE441 ready at: {}", path.display());
//! ```
//!
//! ## References
//!
//! - The Astropy Collaboration (2022). "The Astropy Project: Sustaining
//!   and Growing a Community-oriented Open-Source Project and the Latest
//!   Major Release (v5.0) of the Core Package". *The Astrophysical
//!   Journal* **935**, 167. doi:10.3847/1538-4357/ac7c74.
//! - Acton, C. H., Bachman, N., Semenov, B., Wright, E. (2018). "A look
//!   towards the future in the handling of space science mission
//!   geometry". *Planetary and Space Science* **150**, 9–12.
//!   doi:10.1016/j.pss.2017.02.013.

#[cfg(feature = "runtime-data")]
pub mod runtime;

/// Identifies a scientific dataset known to `siderust`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DatasetId {
    /// JPL DE440 planetary ephemeris (~120 MB BSP).
    De440,
    /// JPL DE441 planetary ephemeris, part 2 (~1.65 GB BSP).
    De441,
    /// VSOP87A planetary theory tables (compiled-in).
    Vsop87a,
    /// VSOP87E planetary theory tables, barycentric (compiled-in).
    Vsop87e,
    /// ELP2000-82B lunar theory tables (compiled-in).
    Elp2000,
    /// IERS Earth Orientation Parameters (owned by the `tempoch` crate).
    IersEop,
}

impl DatasetId {
    /// Short string identifier used in filenames and log messages.
    pub const fn as_str(&self) -> &'static str {
        match self {
            DatasetId::De440 => "de440",
            DatasetId::De441 => "de441",
            DatasetId::Vsop87a => "vsop87a",
            DatasetId::Vsop87e => "vsop87e",
            DatasetId::Elp2000 => "elp2000",
            DatasetId::IersEop => "iers_eop",
        }
    }
}

impl std::fmt::Display for DatasetId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

// ── Acquisition model ────────────────────────────────────────────────────────

/// Metadata for a dataset that is downloaded at runtime.
#[derive(Debug, Clone)]
pub struct RuntimeDownloadMeta {
    /// Primary download URL.
    pub url: &'static str,
    /// Filename on disk inside the cache directory.
    pub filename: &'static str,
    /// Expected SHA-256 hex digest (empty string = skip verification).
    pub sha256: &'static str,
    /// Minimum plausible file size in bytes.
    pub min_size: u64,
    /// Human-readable size hint for progress messages.
    pub size_hint: &'static str,
}

/// How a dataset is acquired.
#[derive(Debug, Clone)]
pub enum Acquisition {
    /// Data is compiled directly into the binary; no external fetch needed.
    Embedded,
    /// Data must be downloaded to the local cache at runtime.
    RuntimeDownload(RuntimeDownloadMeta),
    /// Data is owned and exposed by another crate.
    ExternalProvider {
        /// Crate that owns the data.
        owner_crate: &'static str,
        /// Fully-qualified access path inside that crate.
        access_path: &'static str,
    },
}

/// Static metadata describing a single dataset.
#[derive(Debug, Clone)]
pub struct DatasetMeta {
    /// Dataset identifier.
    pub id: DatasetId,
    /// Human-readable name.
    pub name: &'static str,
    /// How this dataset is acquired.
    pub acquisition: Acquisition,
}

// ── Catalog ──────────────────────────────────────────────────────────────────

/// All datasets known to `siderust`, indexed by [`DatasetId`] ordinal.
///
/// Every [`DatasetId`] variant has exactly one entry here; use [`lookup`] for
/// safe access by id.
pub static DATASETS: &[DatasetMeta] = &[
    // Ordinal 0 — De440
    DatasetMeta {
        id: DatasetId::De440,
        name: "JPL DE440",
        acquisition: Acquisition::RuntimeDownload(RuntimeDownloadMeta {
            url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
            filename: "de440.bsp",
            sha256: "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
            min_size: 100_000_000,
            size_hint: "~120 MB",
        }),
    },
    // Ordinal 1 — De441
    DatasetMeta {
        id: DatasetId::De441,
        name: "JPL DE441 (part 2)",
        acquisition: Acquisition::RuntimeDownload(RuntimeDownloadMeta {
            url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
            filename: "de441_part-2.bsp",
            sha256: "3abb17dae2d78dd34880377544aacb54892104a0d4462b322cb9f4454d4887f6",
            min_size: 1_500_000_000,
            size_hint: "~1.65 GB",
        }),
    },
    // Ordinal 2 — Vsop87a
    DatasetMeta {
        id: DatasetId::Vsop87a,
        name: "VSOP87A planetary theory",
        acquisition: Acquisition::Embedded,
    },
    // Ordinal 3 — Vsop87e
    DatasetMeta {
        id: DatasetId::Vsop87e,
        name: "VSOP87E planetary theory (barycentric)",
        acquisition: Acquisition::Embedded,
    },
    // Ordinal 4 — Elp2000
    DatasetMeta {
        id: DatasetId::Elp2000,
        name: "ELP2000-82B lunar theory",
        acquisition: Acquisition::Embedded,
    },
    // Ordinal 5 — IersEop
    DatasetMeta {
        id: DatasetId::IersEop,
        name: "IERS EOP finals2000A",
        acquisition: Acquisition::ExternalProvider {
            owner_crate: "tempoch",
            access_path: "siderust::astro::eop::IersEop",
        },
    },
];

/// Look up dataset metadata by id.
///
/// The exhaustive `match` ensures a compile error if a new [`DatasetId`]
/// variant is added without a corresponding [`DATASETS`] entry.
pub fn lookup(id: DatasetId) -> &'static DatasetMeta {
    let idx = match id {
        DatasetId::De440 => 0,
        DatasetId::De441 => 1,
        DatasetId::Vsop87a => 2,
        DatasetId::Vsop87e => 3,
        DatasetId::Elp2000 => 4,
        DatasetId::IersEop => 5,
    };
    &DATASETS[idx]
}

// ── Error type ───────────────────────────────────────────────────────────────

/// Errors that can occur during dataset operations.
#[derive(Debug)]
pub enum DatasetError {
    /// I/O error (file system, permissions, etc.).
    Io(std::io::Error),
    /// HTTP download error (requires `runtime-data` feature).
    #[cfg(feature = "runtime-data")]
    Download(String),
    /// Data integrity check failed (checksum mismatch, wrong size, etc.).
    Integrity(String),
    /// The requested dataset cannot be managed by [`crate::data::runtime::DatasetManager`]
    /// because its acquisition kind is not [`Acquisition::RuntimeDownload`].
    NotDownloadable(DatasetId),
    /// SPICE format parse error forwarded from [`crate::formats::spice`].
    Spice(crate::formats::spice::SpiceError),
}

impl std::fmt::Display for DatasetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DatasetError::Io(e) => write!(f, "I/O error: {}", e),
            #[cfg(feature = "runtime-data")]
            DatasetError::Download(msg) => write!(f, "download error: {}", msg),
            DatasetError::Integrity(msg) => write!(f, "integrity check failed: {}", msg),
            DatasetError::NotDownloadable(id) => {
                write!(f, "dataset '{}' is not runtime-downloadable", id)
            }
            DatasetError::Spice(e) => write!(f, "SPICE parse error: {}", e),
        }
    }
}

impl std::error::Error for DatasetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DatasetError::Io(e) => Some(e),
            DatasetError::Spice(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for DatasetError {
    fn from(e: std::io::Error) -> Self {
        DatasetError::Io(e)
    }
}

impl From<crate::formats::spice::SpiceError> for DatasetError {
    fn from(e: crate::formats::spice::SpiceError) -> Self {
        DatasetError::Spice(e)
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── lookup ────────────────────────────────────────────────────────────────

    #[test]
    fn lookup_returns_correct_entry_for_every_variant() {
        assert_eq!(lookup(DatasetId::De440).id, DatasetId::De440);
        assert_eq!(lookup(DatasetId::De441).id, DatasetId::De441);
        assert_eq!(lookup(DatasetId::Vsop87a).id, DatasetId::Vsop87a);
        assert_eq!(lookup(DatasetId::Vsop87e).id, DatasetId::Vsop87e);
        assert_eq!(lookup(DatasetId::Elp2000).id, DatasetId::Elp2000);
        assert_eq!(lookup(DatasetId::IersEop).id, DatasetId::IersEop);
    }

    #[test]
    fn every_dataset_id_present_in_datasets() {
        let all_ids = [
            DatasetId::De440,
            DatasetId::De441,
            DatasetId::Vsop87a,
            DatasetId::Vsop87e,
            DatasetId::Elp2000,
            DatasetId::IersEop,
        ];
        for id in all_ids {
            assert!(
                DATASETS.iter().any(|d| d.id == id),
                "DatasetId::{:?} missing from DATASETS",
                id
            );
        }
    }

    // ── acquisition classification ────────────────────────────────────────────

    #[test]
    fn de440_and_de441_are_runtime_download() {
        for id in [DatasetId::De440, DatasetId::De441] {
            assert!(
                matches!(lookup(id).acquisition, Acquisition::RuntimeDownload(_)),
                "{:?} should be RuntimeDownload",
                id
            );
        }
    }

    #[test]
    fn vsop87a_vsop87e_elp2000_are_embedded() {
        for id in [DatasetId::Vsop87a, DatasetId::Vsop87e, DatasetId::Elp2000] {
            assert!(
                matches!(lookup(id).acquisition, Acquisition::Embedded),
                "{:?} should be Embedded",
                id
            );
        }
    }

    #[test]
    fn iers_eop_is_external_provider() {
        let meta = lookup(DatasetId::IersEop);
        match &meta.acquisition {
            Acquisition::ExternalProvider {
                owner_crate,
                access_path,
            } => {
                assert_eq!(*owner_crate, "tempoch");
                assert!(!access_path.is_empty());
            }
            other => panic!("Expected ExternalProvider, got {:?}", other),
        }
    }

    // ── DatasetId helpers ─────────────────────────────────────────────────────

    #[test]
    fn dataset_id_as_str() {
        assert_eq!(DatasetId::De440.as_str(), "de440");
        assert_eq!(DatasetId::De441.as_str(), "de441");
        assert_eq!(DatasetId::Vsop87a.as_str(), "vsop87a");
        assert_eq!(DatasetId::Vsop87e.as_str(), "vsop87e");
        assert_eq!(DatasetId::Elp2000.as_str(), "elp2000");
        assert_eq!(DatasetId::IersEop.as_str(), "iers_eop");
    }

    #[test]
    fn dataset_id_display() {
        assert_eq!(format!("{}", DatasetId::De440), "de440");
        assert_eq!(format!("{}", DatasetId::De441), "de441");
        assert_eq!(format!("{}", DatasetId::IersEop), "iers_eop");
    }

    #[test]
    fn dataset_id_eq_and_clone() {
        let id = DatasetId::De440;
        let id2 = id;
        assert_eq!(id, id2);
        assert_ne!(id, DatasetId::De441);
    }

    // ── DatasetError ──────────────────────────────────────────────────────────

    #[test]
    fn dataset_error_display_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let e = DatasetError::Io(io_err);
        assert!(format!("{}", e).contains("I/O error"));
    }

    #[test]
    fn dataset_error_display_integrity() {
        let e = DatasetError::Integrity("checksum mismatch".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("integrity check failed"));
        assert!(msg.contains("checksum mismatch"));
    }

    #[test]
    fn dataset_error_display_not_downloadable() {
        let e = DatasetError::NotDownloadable(DatasetId::Vsop87a);
        let msg = format!("{}", e);
        assert!(msg.contains("vsop87a"), "msg: {msg}");
        assert!(msg.contains("not runtime-downloadable"), "msg: {msg}");
    }

    #[test]
    fn dataset_error_display_spice() {
        let e = DatasetError::Spice(crate::formats::spice::SpiceError::FormatParse(
            "bad header".to_string(),
        ));
        let msg = format!("{}", e);
        assert!(msg.contains("SPICE"), "msg: {msg}");
        assert!(msg.contains("bad header"), "msg: {msg}");
    }

    #[test]
    fn dataset_error_from_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "denied");
        let e: DatasetError = io_err.into();
        assert!(matches!(e, DatasetError::Io(_)));
    }

    #[test]
    fn dataset_error_from_spice() {
        let spice_err = crate::formats::spice::SpiceError::FormatParse("test".to_string());
        let e: DatasetError = spice_err.into();
        assert!(matches!(e, DatasetError::Spice(_)));
    }

    #[test]
    fn dataset_error_source_io_is_some() {
        use std::error::Error;
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "not found");
        let e = DatasetError::Io(io_err);
        assert!(e.source().is_some());
    }

    #[test]
    fn dataset_error_source_non_io_is_none() {
        use std::error::Error;
        let e = DatasetError::Integrity("x".to_string());
        assert!(e.source().is_none());
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn dataset_error_display_download() {
        let e = DatasetError::Download("timeout".to_string());
        let msg = format!("{}", e);
        assert!(msg.contains("download error"));
        assert!(msg.contains("timeout"));
    }

    // ── RuntimeDownloadMeta fields ────────────────────────────────────────────

    #[test]
    fn runtime_download_meta_fields_non_empty_for_de440() {
        let meta = lookup(DatasetId::De440);
        if let Acquisition::RuntimeDownload(rdm) = &meta.acquisition {
            assert!(!rdm.url.is_empty());
            assert!(!rdm.filename.is_empty());
            assert!(rdm.min_size > 0);
            assert!(!rdm.size_hint.is_empty());
        } else {
            panic!("De440 must be RuntimeDownload");
        }
    }
}
