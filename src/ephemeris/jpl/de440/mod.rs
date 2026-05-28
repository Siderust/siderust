// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![cfg_attr(siderust_mock_de440, allow(dead_code))]

//! # DE440 Ephemeris
//!
//! ## Scientific scope
//!
//! DE440 is the JPL planetary and lunar ephemeris released in 2021.  It is
//! fitted to modern observational data (1950–2050 range of high-accuracy
//! measurements) including lunar laser ranging (LLR) and Mars-lander tracking.
//! Time span: JD 2287184.5 – 2688976.5 (1549-Dec-31 to 2650-Jan-25).
//! Coordinate frame: ICRF (km, TDB).
//!
//! ## Technical scope
//!
//! Thin configuration layer over the shared JPL DE4xx infrastructure.
//! Provides [`De440Data`], the marker type that selects DE440 coefficients
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
use crate::ephemeris::jpl::jpl_data::de440 as data;

/// Marker type selecting DE440 embedded coefficient data.
pub(crate) struct De440Data;

impl DeData for De440Data {
    const SUN: SegmentDescriptor = data::SUN;
    const EMB: SegmentDescriptor = data::EMB;
    const MOON: SegmentDescriptor = data::MOON;
}
