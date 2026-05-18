// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Binary file-format parsers
//!
//! Low-level parsers for standard astronomical data formats.
//! These modules operate on raw byte slices without knowledge of
//! the dataset catalog or acquisition machinery.
//!
//! For the dataset catalog (what datasets exist and how to acquire them) see
//! [`siderust::datasets`](crate::datasets). The `formats` modules sit *below*
//! the catalog: they are called by the runtime back-end after a file has
//! been located on disk.

pub mod spice;
