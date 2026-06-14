// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # In-memory large-star catalog
//!
//! [`LargeStarCatalog`] is the in-memory container with cone search
//! (`cone_search`), magnitude/quality filters, and an epoch-aware
//! query that propagates positions via proper motion.

use crate::qtty::Radians;

use super::record::{inside_cone, passes_filter, CatalogFilter, CatalogRecord};

/// In-memory large-star catalog.
#[derive(Debug, Default, Clone)]
pub struct LargeStarCatalog {
    records: Vec<CatalogRecord>,
}

impl LargeStarCatalog {
    /// Empty catalog.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of records held.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// True iff the catalog has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Append a single record.
    pub fn push(&mut self, record: CatalogRecord) {
        self.records.push(record);
    }

    /// Extend with an iterator of records.
    pub fn extend(&mut self, iter: impl IntoIterator<Item = CatalogRecord>) {
        self.records.extend(iter);
    }

    /// Borrow the underlying record slice.
    pub fn records(&self) -> &[CatalogRecord] {
        &self.records
    }

    /// Cone search around the supplied ICRS direction. Records are evaluated at
    /// their catalog epoch (no proper-motion shift); use
    /// [`Self::cone_search_at_epoch`] for epoch-aware queries.
    pub fn cone_search(
        &self,
        center_ra: Radians,
        center_dec: Radians,
        radius: Radians,
        filter: CatalogFilter,
    ) -> Vec<&CatalogRecord> {
        self.records
            .iter()
            .filter(|r| {
                inside_cone(r.ra, r.dec, center_ra, center_dec, radius) && passes_filter(r, &filter)
            })
            .collect()
    }

    /// Cone search at a specified epoch. Each record is propagated to
    /// `target_epoch_jyr` via proper motion before the angular-distance test
    /// runs. Records without proper motion are tested at their catalog epoch.
    pub fn cone_search_at_epoch(
        &self,
        center_ra: Radians,
        center_dec: Radians,
        radius: Radians,
        target_epoch_jyr: f64,
        filter: CatalogFilter,
    ) -> Vec<&CatalogRecord> {
        self.records
            .iter()
            .filter(|r| {
                let (ra, dec) = r.propagate_to(target_epoch_jyr);
                inside_cone(ra, dec, center_ra, center_dec, radius) && passes_filter(r, &filter)
            })
            .collect()
    }

    /// Chunked cone search: splits the underlying record slice into
    /// `chunk_size`-row windows and yields matches per chunk. The concatenated
    /// result is identical to [`Self::cone_search`].
    pub fn cone_search_chunked(
        &self,
        center_ra: Radians,
        center_dec: Radians,
        radius: Radians,
        filter: CatalogFilter,
        chunk_size: usize,
    ) -> impl Iterator<Item = Vec<&CatalogRecord>> {
        let chunks = self.records.chunks(chunk_size.max(1));
        chunks.map(move |chunk| {
            chunk
                .iter()
                .filter(|r| {
                    inside_cone(r.ra, r.dec, center_ra, center_dec, radius)
                        && passes_filter(r, &filter)
                })
                .collect()
        })
    }
}
