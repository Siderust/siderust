// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU 2000A nutation model marker.
//!
//! This selects the full MHB2000 nutation series with 1365 trigonometric
//! terms and no IAU 2006 compatibility correction.

use super::{Model, NutationModelId, NutationTag};

/// Full IAU 2000A nutation model.
///
/// Use this when you want the uncorrected IAU 2000A series rather than the
/// default IAU 2006A-compatible variant.
pub type Iau2000A = Model<Tag>;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Tag;

impl NutationTag for Tag {
    const ID: NutationModelId = NutationModelId::Iau2000A;
}
