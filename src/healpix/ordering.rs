// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

/// HEALPix ordering scheme.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum HealpixOrdering {
    /// HEALPix RING ordering.
    Ring,
    /// HEALPix NESTED ordering.
    Nested,
}
