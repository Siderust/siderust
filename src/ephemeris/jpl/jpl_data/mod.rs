// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! JPL planetary ephemeris data sub-modules.
//!
//! Each child module exposes the `SegmentDescriptor` constants (`SUN`, `EMB`,
//! `MOON`) that the corresponding ephemeris backend consumes at compile time.

#[cfg(feature = "de440")]
pub(crate) mod de440;
#[cfg(feature = "de441")]
pub(crate) mod de441;
