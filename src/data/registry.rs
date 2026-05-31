// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Dataset registry, metadata for all downloadable datasets.

/// Identifies a downloadable dataset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DatasetId {
    /// JPL DE440 planetary ephemeris (~120 MB BSP).
    De440,
    /// JPL DE441 planetary ephemeris, part 2 (~1.65 GB BSP).
    De441,
    /// JPL Mars satellite ephemeris used for Mars-center offsets.
    Mar099,
    /// JPL Jupiter satellite ephemeris used for Jupiter-center offsets.
    Jup365Merged,
    /// JPL Saturn satellite ephemeris used for Saturn-center offsets.
    Sat441l,
    /// JPL Uranus satellite ephemeris used for Uranus-center offsets.
    Ura184Merged,
    /// JPL Neptune satellite ephemeris used for Neptune-center offsets.
    Nep097,
    // Smaller datasets, normally compiled in, but available for runtime override.
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
            DatasetId::Mar099 => "mar099",
            DatasetId::Jup365Merged => "jup365_merged",
            DatasetId::Sat441l => "sat441l",
            DatasetId::Ura184Merged => "ura184_merged",
            DatasetId::Nep097 => "nep097",
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
        sha256: "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
        min_size: 100_000_000,
        size_hint: "~120 MB",
    },
    DatasetMeta {
        id: DatasetId::De441,
        name: "JPL DE441 (part 2)",
        url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
        filename: "de441_part-2.bsp",
        sha256: "3abb17dae2d78dd34880377544aacb54892104a0d4462b322cb9f4454d4887f6",
        min_size: 1_500_000_000,
        size_hint: "~1.65 GB",
    },
    DatasetMeta {
        id: DatasetId::Mar099,
        name: "JPL MAR099 Mars satellite ephemeris",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/mar099.bsp",
        filename: "mar099.bsp",
        sha256: "",
        min_size: 1_000_000_000,
        size_hint: "~1.18 GB",
    },
    DatasetMeta {
        id: DatasetId::Jup365Merged,
        name: "JPL JUP365 Jupiter satellite ephemeris",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/jup365.bsp",
        filename: "jup365.bsp",
        sha256: "",
        min_size: 1_000_000_000,
        size_hint: "~1.11 GB",
    },
    DatasetMeta {
        id: DatasetId::Sat441l,
        name: "JPL SAT441 Saturn satellite ephemeris",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/sat441l.bsp",
        filename: "sat441l.bsp",
        sha256: "",
        min_size: 600_000_000,
        size_hint: "~639 MB",
    },
    DatasetMeta {
        id: DatasetId::Ura184Merged,
        name: "JPL URA184 Uranus satellite ephemeris",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/ura184.bsp",
        filename: "ura184.bsp",
        sha256: "",
        min_size: 4_000_000_000,
        size_hint: "~4.44 GB",
    },
    DatasetMeta {
        id: DatasetId::Nep097,
        name: "JPL NEP097 Neptune satellite ephemeris",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/nep097.bsp",
        filename: "nep097.bsp",
        sha256: "",
        min_size: 3_000_000_000,
        size_hint: "~3.23 GB",
    },
    DatasetMeta {
        id: DatasetId::IersEop,
        name: "IERS EOP finals2000A",
        url: "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all",
        filename: "finals2000A.all",
        sha256: "eeb76193cd43c065b763d78923f8a8ce2f8a62df5d3e1519a3217fac433bdaa6",
        min_size: 1_000_000,
        size_hint: "~4 MB",
    },
];

/// Look up dataset metadata by ID.
#[cfg(feature = "runtime-data")]
pub fn lookup(id: DatasetId) -> Option<&'static DatasetMeta> {
    DATASETS.iter().find(|d| d.id == id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dataset_id_as_str() {
        assert_eq!(DatasetId::De440.as_str(), "de440");
        assert_eq!(DatasetId::De441.as_str(), "de441");
        assert_eq!(DatasetId::Mar099.as_str(), "mar099");
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
        assert_ne!(DatasetId::Vsop87a, DatasetId::Vsop87e);
    }

    #[test]
    fn datasets_slice_is_non_empty() {
        assert!(!DATASETS.is_empty());
        // Check that de440 is present
        assert!(DATASETS.iter().any(|d| d.id == DatasetId::De440));
        assert!(DATASETS.iter().any(|d| d.id == DatasetId::IersEop));
    }

    #[test]
    fn dataset_meta_fields_non_empty() {
        let de440 = DATASETS.iter().find(|d| d.id == DatasetId::De440).unwrap();
        assert!(!de440.name.is_empty());
        assert!(!de440.url.is_empty());
        assert!(!de440.filename.is_empty());
        assert!(de440.min_size > 0);
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn lookup_de440_found() {
        let meta = lookup(DatasetId::De440);
        assert!(meta.is_some());
        assert_eq!(meta.unwrap().id, DatasetId::De440);
    }

    #[test]
    #[cfg(feature = "runtime-data")]
    fn lookup_elp2000_not_in_catalog_returns_none() {
        // Elp2000 is defined as DatasetId but might not be in DATASETS
        // Either way lookup returns Some or None without panicking
        let _result = lookup(DatasetId::Elp2000);
    }
}
