// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Event searches
//!
//! ## Scientific scope
//! Altitude/azimuth event searches, solar/lunar/stellar rise-transit-set,
//! and period-finding for astronomical observation planning.
//!
//! ## Technical scope
//! Entry points: [`altitude`], [`azimuth`], [`solar`], [`lunar`], [`stellar`].
//! Internal coordinate support in [`horizontal`].

pub mod altitude;
pub mod azimuth;
pub mod horizontal;
pub mod lunar;
pub mod solar;
pub mod stellar;
