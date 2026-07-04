// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Gaia DR3 quality metadata carried with raw and validated source records.

/// Minimal Gaia DR3 quality flags preserved by the typed catalogue layer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GaiaDr3QualityFlags {
    /// Whether the row passed caller-selected quality filtering before ingestion.
    pub quality_ok: bool,
    /// Gaia `duplicated_source`, when present.
    pub duplicated_source: Option<bool>,
}

impl Default for GaiaDr3QualityFlags {
    fn default() -> Self {
        Self {
            quality_ok: true,
            duplicated_source: None,
        }
    }
}
