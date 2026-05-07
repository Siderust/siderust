// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2000A Nutation Model Marker
//!
//! Type-level marker selecting the full uncorrected IAU 2000A nutation
//! series in the transform pipeline.
//!
//! ## Scientific scope
//!
//! The IAU 2000A nutation series (MHB2000) consists of 1365 trigonometric
//! terms — 678 luni-solar plus 687 planetary — and delivers sub-microarcsecond
//! formal accuracy for the nutation angles `Δψ` and `Δε`. This marker
//! selects 2000A *without* the Wallace-Capitaine P03/J₂ correction that
//! makes it strictly compatible with IAU 2006 precession; use [`Iau2006A`]
//! for the corrected variant.
//!
//! ## Technical scope
//!
//! [`Iau2000A`] is an alias for `Model<Tag>` where `Tag` is a sealed,
//! zero-sized type implementing [`NutationTag`] with
//! `ID = NutationModelId::Iau2000A`. The actual numerical evaluation lives
//! in the shared `nut00a` engine; this file only fixes the compile-time
//! identity used for trait dispatch.
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.6
//! * Mathews, Herring & Buffett (2002), *J. Geophys. Res.* 107, B4
//! * SOFA routine `iauNut00a`

use super::{Model, NutationModelId, NutationTag};

/// Full IAU 2000A nutation model.
///
/// Use this when you want the uncorrected IAU 2000A series rather than the
/// default IAU 2006A-compatible variant.
pub type Iau2000A = Model<Tag>;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Tag;

impl super::private::Sealed for Tag {}

impl NutationTag for Tag {
    const ID: NutationModelId = NutationModelId::Iau2000A;
}
