// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # DE440 Ephemeris Module
//!
//! Thin configuration layer over the shared JPL DE4xx infrastructure.
//! Provides [`De440Data`] — the marker type that selects DE440 coefficients
//! for use with [`DeEphemeris`](crate::calculus::jpl::DeEphemeris).
//!
//! All evaluation logic, body-chain arithmetic, and frame/unit conversions
//! are handled by [`calculus::jpl`](crate::calculus::jpl).
//!
//! ## Reference Frame & Units
//!
//! DE440 data is in ICRF (≈ J2000 equatorial), km and km/s, on TDB.
//! The public API converts to ecliptic AU via the `Ephemeris` trait.

pub mod data;

use super::{eval::SegmentDescriptor, DeData};

/// Marker type selecting DE440 embedded coefficient data.
pub struct De440Data;

impl DeData for De440Data {
    const SUN: SegmentDescriptor = data::SUN;
    const EMB: SegmentDescriptor = data::EMB;
    const MOON: SegmentDescriptor = data::MOON;
}
