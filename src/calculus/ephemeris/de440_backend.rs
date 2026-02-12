// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE440 ephemeris backend.
//!
//! Thin alias over the shared JPL DE4xx backend implementation.

use crate::calculus::jpl::{de440::De440Data, DeEphemeris};

/// JPL DE440 ephemeris backend.
pub type De440Ephemeris = DeEphemeris<De440Data>;
