// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE441 ephemeris backend.
//!
//! Thin alias over the shared JPL DE4xx backend implementation.

#[cfg(not(siderust_mock_de441))]
use crate::calculus::jpl::{de441::De441Data, DeEphemeris};

/// JPL DE441 ephemeris backend.
///
/// When built with `SIDERUST_JPL_STUB=de441` (or `all`), this type aliases to
/// [`Vsop87Ephemeris`](super::Vsop87Ephemeris) as a mock backend so tests can
/// run without downloading the DE441 BSP.
#[cfg(not(siderust_mock_de441))]
pub type De441Ephemeris = DeEphemeris<De441Data>;

/// Mock DE441 backend used for no-download builds/tests.
#[cfg(siderust_mock_de441)]
pub type De441Ephemeris = super::Vsop87Ephemeris;
