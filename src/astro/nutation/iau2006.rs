// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2006 precession-only model marker.
//!
//! This disables nutation offsets entirely while keeping the IAU 2006 mean
//! obliquity, which is useful for comparisons against precession-only chains.

use super::{Model, NutationModelId, NutationTag};

/// IAU 2006 precession-only profile.
///
/// The returned nutation angles are zero, but the IAU 2006 mean obliquity is
/// still used so frame paths remain internally consistent.
pub type Iau2006 = Model<Tag>;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Tag;

impl super::private::Sealed for Tag {}

impl NutationTag for Tag {
    const ID: NutationModelId = NutationModelId::Iau2006;
}
