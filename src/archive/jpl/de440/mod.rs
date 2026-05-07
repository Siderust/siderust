// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded DE440 coefficient data.
//!
//! The build script generates `de440_data.rs` which defines sub-modules
//! `sun`, `emb`, and `moon`, plus pre-built `SUN`, `EMB`, `MOON`
//! [`SegmentDescriptor`] constants.

use crate::calculus::jpl::eval::SegmentDescriptor;
#[allow(unused_imports)]
use crate::qtty;

include!(concat!(env!("OUT_DIR"), "/de440_data.rs"));
