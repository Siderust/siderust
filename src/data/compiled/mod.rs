// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Crate-private compiled data tables and generated coefficient archives.

pub(crate) mod jpl;
#[cfg(feature = "lagrange-centers")]
pub(crate) mod lagrange;
pub(crate) mod nut00a_tables;
pub(crate) mod pluto_tables;
