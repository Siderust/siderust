// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![cfg_attr(siderust_mock_de440, allow(dead_code))]
#![allow(unreachable_pub)]

//! Embedded DE440 coefficient data.
//!
//! The build script generates `de440_data.rs` which defines sub-modules
//! `sun`, `emb`, and `moon`, plus pre-built `SUN`, `EMB`, `MOON`
//! [`SegmentDescriptor`] constants.

use crate::ephemeris::jpl::eval::SegmentDescriptor;
#[allow(unused_imports)]
use crate::qtty;

include!(concat!(env!("OUT_DIR"), "/de440_data.rs"));
