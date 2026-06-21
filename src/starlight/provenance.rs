// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

/// Provenance metadata for generated stellar maps.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StellarMapProvenance {
    /// Dataset name of the generated map.
    pub dataset_name: String,
    /// Dataset version string.
    pub version: String,
    /// UTC generation date or timestamp.
    pub generation_date_utc: String,
    /// Source catalogue name or path.
    pub source_catalogue: String,
    /// Optional source catalogue release identifier.
    pub source_catalogue_release: Option<String>,
    /// Optional source catalogue license.
    pub source_catalogue_license: Option<String>,
    /// Optional source catalogue checksum.
    pub source_catalogue_checksum: Option<String>,
    /// Optional documented magnitude limit.
    pub magnitude_limit: Option<String>,
    /// Photometric band definition.
    pub band_definition: String,
    /// Photometry model identifier.
    pub photometry_model: String,
    /// Optional smoothing description.
    pub smoothing: Option<String>,
    /// Generator identifier.
    pub generator: String,
}
