// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Binary file-format parsers
//!
//! Low-level parsers for standard astronomical data formats.
//! These modules operate on raw byte slices without knowledge of
//! the dataset catalog or acquisition machinery.

pub mod spice;
