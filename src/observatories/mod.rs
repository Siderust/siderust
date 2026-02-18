// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory Catalog Module
//!
//! This module provides constant definitions for well-known observatories,
//! stored as [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic) values
//! (geodetic longitude, latitude, and ellipsoidal height above the WGS84
//! ellipsoid associated with the [`ECEF`](crate::coordinates::frames::ECEF)
//! frame).
//!
//! ## Geodetic vs Spherical
//!
//! Observatory constants are intentionally stored as ellipsoidal positions,
//! **not** as spherical `Position` objects.  A geodetic position
//! (lon/lat/height-above-ellipsoid) cannot be round-tripped correctly through a
//! spherical position's `.to_cartesian()` because the spherical `distance` field
//! represents a **radial distance**, not an ellipsoidal height.  The correct
//! conversion from geodetic to ECEF requires an ellipsoid-aware formula
//! (WGS84 Bowring formula), available via
//! [`to_cartesian`](affn::ellipsoidal::Position::to_cartesian).
//!
//! ## Usage
//!
//! To obtain the WGS84 ECEF Cartesian position of an observatory:
//!
//! ```rust
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::ECEF;
//! use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
//! use qtty::Meter;
//!
//! let ecef = ROQUE_DE_LOS_MUCHACHOS.to_cartesian::<Meter>();
//! ```

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use qtty::{Degrees, Meters};

/// Roque de los Muchachos Observatory (La Palma, Canary Islands, Spain)
/// - Longitude: −17.8925°
/// - Latitude:  +28.7543°
/// - Altitude:  2396 m (WGS84 ellipsoidal height)
pub const ROQUE_DE_LOS_MUCHACHOS: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-17.8925),
    Degrees::new(28.7543),
    Meters::new(2396.0),
);

/// El Paranal (Chile)
/// - Longitude: −70.4043°
/// - Latitude:  −24.6272°
/// - Altitude:  2635 m (WGS84 ellipsoidal height)
pub const EL_PARANAL: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-70.4043),
    Degrees::new(-24.6272),
    Meters::new(2635.0),
);

/// Mauna Kea Observatory (Hawaiʻi, USA)
/// - Longitude: −155.4681°
/// - Latitude:  +19.8207°
/// - Altitude:  4207 m (WGS84 ellipsoidal height)
pub const MAUNA_KEA: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-155.4681),
    Degrees::new(19.8207),
    Meters::new(4207.0),
);

/// La Silla Observatory (ESO, Chile)
/// - Longitude: −70.7346°
/// - Latitude:  −29.2584°
/// - Altitude:  2400 m (WGS84 ellipsoidal height)
pub const LA_SILLA_OBSERVATORY: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-70.7346),
    Degrees::new(-29.2584),
    Meters::new(2400.0),
);
