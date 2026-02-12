// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE441 ephemeris backend.
//!
//! Thin alias over the shared JPL DE4xx backend implementation.

use crate::calculus::jpl::{de441::De441Data, DeEphemeris};

/// JPL DE441 ephemeris backend.
pub type De441Ephemeris = DeEphemeris<De441Data>;
