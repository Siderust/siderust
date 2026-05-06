// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006 Precession-only Nutation Marker
//!
//! Type-level marker selecting a "precession-only" profile: the IAU 2006
//! mean obliquity is retained but the nutation angles are forced to zero.
//!
//! ## Scientific scope
//!
//! For diagnostic comparisons it is often useful to evaluate the
//! orientation chain with precession only and no nutation, isolating the
//! secular drift of the equator from the periodic wobble. Forcing
//! `Δψ = Δε = 0` while keeping the IAU 2006 mean obliquity `ε_A` ensures
//! that the rest of the chain (CIO locator, ERA, polar motion) remains
//! internally consistent and only the nutation contribution is suppressed.
//!
//! ## Technical scope
//!
//! [`Iau2006`] is an alias for `Model<Tag>` whose `Tag` implements the
//! sealed [`NutationTag`] trait with
//! `ID = NutationModelId::Iau2006`. Trait dispatch routes to a degenerate
//! evaluator that returns zero nutation while still reporting the
//! correctly evaluated mean obliquity.
//!
//! ## References
//!
//! * IAU 2006 Resolution B1
//! * IERS Conventions (2010), §5.6

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
