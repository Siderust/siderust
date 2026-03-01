// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE440 ephemeris backend.
//!
//! Thin alias over the shared JPL DE4xx backend implementation.

#[cfg(not(siderust_mock_de440))]
use crate::calculus::jpl::{de440::De440Data, DeEphemeris};

/// JPL DE440 ephemeris backend.
///
/// When built with `SIDERUST_JPL_STUB=de440` (or `all`), this type aliases to
/// [`Vsop87Ephemeris`](super::Vsop87Ephemeris) as a mock backend so tests can
/// run without downloading the DE440 BSP.
#[cfg(not(siderust_mock_de440))]
pub type De440Ephemeris = DeEphemeris<De440Data>;

/// Mock DE440 backend used for no-download builds/tests.
#[cfg(siderust_mock_de440)]
pub type De440Ephemeris = super::Vsop87Ephemeris;
