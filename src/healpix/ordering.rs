// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! HEALPix pixel ordering schemes.
//!
//! RING ordering stores pixels by iso-latitude rings. NESTED ordering stores
//! pixels according to a hierarchical quadtree-like subdivision and requires a
//! power-of-two `nside`.

/// HEALPix ordering scheme.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum HealpixOrdering {
    /// HEALPix RING ordering by iso-latitude rings.
    Ring,
    /// HEALPix NESTED hierarchical ordering.
    Nested,
}
