// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Dataset registry — metadata for all downloadable datasets.

/// Identifies a downloadable dataset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DatasetId {
    /// JPL DE440 planetary ephemeris (~120 MB BSP).
    De440,
    /// JPL DE441 planetary ephemeris, part 2 (~1.65 GB BSP).
    De441,
    // Smaller datasets — normally compiled in, but available for runtime override.
    /// VSOP87A planetary theory tables.
    Vsop87a,
    /// VSOP87E planetary theory tables (barycentric).
    Vsop87e,
    /// ELP2000-82B lunar theory tables.
    Elp2000,
    /// IERS Earth Orientation Parameters (finals2000A.all).
    IersEop,
}

impl DatasetId {
    /// Short string identifier for use in filenames and logs.
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

/// Static metadata describing a single downloadable dataset.
#[derive(Debug, Clone)]
pub struct DatasetMeta {
    /// Dataset identifier.
    pub id: DatasetId,
    /// Human-readable name.
    pub name: &'static str,
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

/// All known datasets.
pub static DATASETS: &[DatasetMeta] = &[
    DatasetMeta {
        id: DatasetId::De440,
        name: "JPL DE440",
        url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
        filename: "de440.bsp",
        sha256: "", // TODO: compute and pin
        min_size: 100_000_000,
        size_hint: "~120 MB",
    },
    DatasetMeta {
        id: DatasetId::De441,
        name: "JPL DE441 (part 2)",
        url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
        filename: "de441_part-2.bsp",
        sha256: "", // TODO: compute and pin
        min_size: 1_500_000_000,
        size_hint: "~1.65 GB",
    },
    DatasetMeta {
        id: DatasetId::IersEop,
        name: "IERS EOP finals2000A",
        url: "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all",
        filename: "finals2000A.all",
        sha256: "",
        min_size: 1_000_000,
        size_hint: "~4 MB",
    },
];

/// Look up dataset metadata by ID.
#[cfg(feature = "runtime-data")]
pub fn lookup(id: DatasetId) -> Option<&'static DatasetMeta> {
    DATASETS.iter().find(|d| d.id == id)
}
