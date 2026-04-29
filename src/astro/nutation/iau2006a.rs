// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2006A nutation model marker.
//!
//! This selects the full MHB2000 nutation series together with the
//! Wallace-Capitaine P03/J2 correction so it is compatible with the IAU 2006
//! precession model.

use super::{Model, NutationModelId, NutationTag};

/// IAU 2006-compatible full nutation model.
///
/// This is the default high-precision choice used by the transform pipeline.
pub type Iau2006A = Model<Tag>;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Tag;

impl super::private::Sealed for Tag {}

impl NutationTag for Tag {
    const ID: NutationModelId = NutationModelId::Iau2006A;
}
