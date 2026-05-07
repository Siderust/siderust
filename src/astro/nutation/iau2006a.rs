// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006A Nutation Model Marker
//!
//! Type-level marker selecting the full MHB2000 series with the Wallace-
//! Capitaine P03/J₂ correction so that it is consistent with IAU 2006
//! precession.
//!
//! ## Scientific scope
//!
//! The raw IAU 2000A series was developed against the IAU 1976 precession
//! and exhibits secular residuals when combined with the modernised IAU
//! 2006 (P03) precession. Wallace & Capitaine (2006) derived a small
//! polynomial correction to `Δψ` and `Δε` that removes those residuals and
//! preserves sub-microarcsecond consistency between the two models. The
//! corrected series is the high-precision default for IAU 2006/2000A
//! reductions and is what [`Iau2006A`] selects.
//!
//! ## Technical scope
//!
//! [`Iau2006A`] is an alias for `Model<Tag>`. `Tag` implements the sealed
//! internal marker trait with `ID = NutationModelId::Iau2006A`, dispatching
//! to the IAU 2006A-corrected evaluator inside the shared `nut00a` engine.
//! This is the model used by the default transform pipeline.
//!
//! ## References
//!
//! * Wallace, P. T. & Capitaine, N. (2006), A&A 459, 981
//! * IERS Conventions (2010), §5.5.4
//! * SOFA routine `iauNut06a`

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
