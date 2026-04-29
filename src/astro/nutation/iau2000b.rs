// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2000B nutation model marker.
//!
//! This selects the abridged 77-term nutation series with a fixed planetary
//! correction, matching the public [`nutation_iau2000b`](super::nutation_iau2000b)
//! function.

use super::{Model, NutationModelId, NutationTag};

/// Abridged IAU 2000B nutation model.
///
/// This is the low-cost alternative to the full 2000A/2006A series and
/// remains accurate to roughly the milliarcsecond level.
pub type Iau2000B = Model<Tag>;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Tag;

impl super::private::Sealed for Tag {}

impl NutationTag for Tag {
    const ID: NutationModelId = NutationModelId::Iau2000B;
}
