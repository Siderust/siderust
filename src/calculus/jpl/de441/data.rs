// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded DE441 coefficient data.
//!
//! The build script generates `de441_data.rs` which defines sub-modules
//! `sun`, `emb`, and `moon`, plus pre-built `SUN`, `EMB`, `MOON`
//! [`SegmentDescriptor`] constants.

use super::super::eval::SegmentDescriptor;

include!(concat!(env!("OUT_DIR"), "/de441_data.rs"));
