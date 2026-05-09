// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2000B Nutation Model Marker
//!
//! Type-level marker selecting the abridged IAU 2000B nutation series in
//! the transform pipeline.
//!
//! ## Scientific scope
//!
//! IAU 2000B truncates the full MHB2000 series to **77 luni-solar terms**
//! and replaces the omitted planetary contributions by a fixed correction.
//! It reproduces the full IAU 2000A nutation to better than ≈ 1 mas while
//! costing roughly an order of magnitude less to evaluate, which makes it
//! the standard low-cost choice for ephemeris and pointing pipelines that
//! do not require sub-mas accuracy.
//!
//! ## Technical scope
//!
//! [`Iau2000B`] is an alias for `Model<Tag>` whose `Tag` implements the
//! sealed internal marker trait with
//! `ID = NutationModelId::Iau2000B`, so dispatch in
//! [`AstroContext`](crate::coordinates::transform::context::AstroContext)
//! routes to the abridged evaluator without any runtime branching.
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.6
//! * McCarthy & Luzum (2003), *Celest. Mech.* 85, 37–49
//! * SOFA routine `iauNut00b`

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
