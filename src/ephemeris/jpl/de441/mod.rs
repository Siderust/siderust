// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![cfg_attr(siderust_mock_de441, allow(dead_code))]

//! # DE441 Ephemeris
//!
//! ## Scientific scope
//!
//! DE441 is the JPL extended-range planetary and lunar ephemeris (2021).  It
//! shares the same dynamical model as DE440 but covers a much longer time
//! span at the cost of a slightly looser fit to modern observational data.
//! Time span: JD 625648.5 – 2816912.5 (−13 200 to +17 191, approximately).
//! Coordinate frame: ICRF (km, TDB).
//!
//! ## Technical scope
//!
//! Thin configuration layer over the shared JPL DE4xx infrastructure.
//! Provides [`De441Data`], the marker type that selects DE441 coefficients
//! for use with [`DeEphemeris`](crate::ephemeris::jpl::DeEphemeris).
//!
//! All evaluation logic, body-chain arithmetic, and frame / unit conversions
//! are handled by the internal [`crate::ephemeris::jpl`] implementation.
//!
//! ## References
//!
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides
//!   DE440 and DE441". *The Astronomical Journal* 161, 105.
//!   <https://doi.org/10.3847/1538-3881/abd414>

use super::{eval::SegmentDescriptor, DeData};
use crate::ephemeris::jpl::jpl_data::de441 as data;

/// Marker type selecting DE441 embedded coefficient data.
pub(crate) struct De441Data;

impl DeData for De441Data {
    const SUN: SegmentDescriptor = data::SUN;
    const EMB: SegmentDescriptor = data::EMB;
    const MOON: SegmentDescriptor = data::MOON;
}
