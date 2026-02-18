// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory Catalog Module
//!
//! This module provides constant definitions for well-known observatories,
//! stored as [`GeodeticCoord`] values (geodetic longitude, latitude, and
//! ellipsoidal height above WGS84).
//!
//! ## Geodetic vs Spherical
//!
//! Observatory constants are intentionally stored as [`GeodeticCoord`], **not**
//! as spherical `Position` objects. A geodetic position (lon/lat/height-above-
//! ellipsoid) cannot be round-tripped correctly through a spherical position's
//! `.to_cartesian()` because the spherical `distance` field represents a
//! **radial distance**, not an ellipsoidal height. The correct conversion from
//! geodetic to ECEF requires an ellipsoid-aware formula (such as the WGS84
//! Bowring formula implemented in `ObserverSite::geocentric_itrf()`).
//!
//! ## Usage
//!
//! To obtain the WGS84 ECEF Cartesian position of an observatory:
//!
//! ```rust
//! use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
//! use siderust::coordinates::centers::ObserverSite;
//! use qtty::Meter;
//!
//! let site = ObserverSite::from_geodetic(&ROQUE_DE_LOS_MUCHACHOS);
//! let ecef = site.geocentric_itrf::<Meter>();
//! ```

use affn::geodesy::GeodeticCoord;
use qtty::{Degrees, Meters};

/// Roque de los Muchachos Observatory (La Palma, Canary Islands, Spain)
/// - Longitude: −17.8925°
/// - Latitude:  +28.7543°
/// - Altitude:  2396 m (WGS84 ellipsoidal height)
pub const ROQUE_DE_LOS_MUCHACHOS: GeodeticCoord = GeodeticCoord {
    lon: Degrees::new(-17.8925),
    lat: Degrees::new(28.7543),
    height: Meters::new(2396.0),
};

/// El Paranal (Chile)
/// - Longitude: −70.4043°
/// - Latitude:  −24.6272°
/// - Altitude:  2635 m (WGS84 ellipsoidal height)
pub const EL_PARANAL: GeodeticCoord = GeodeticCoord {
    lon: Degrees::new(-70.4043),
    lat: Degrees::new(-24.6272),
    height: Meters::new(2635.0),
};

/// Mauna Kea Observatory (Hawaiʻi, USA)
/// - Longitude: −155.4681°
/// - Latitude:  +19.8207°
/// - Altitude:  4207 m (WGS84 ellipsoidal height)
pub const MAUNA_KEA: GeodeticCoord = GeodeticCoord {
    lon: Degrees::new(-155.4681),
    lat: Degrees::new(19.8207),
    height: Meters::new(4207.0),
};

/// La Silla Observatory (ESO, Chile)
/// - Longitude: −70.7346°
/// - Latitude:  −29.2584°
/// - Altitude:  2400 m (WGS84 ellipsoidal height)
pub const LA_SILLA_OBSERVATORY: GeodeticCoord = GeodeticCoord {
    lon: Degrees::new(-70.7346),
    lat: Degrees::new(-29.2584),
    height: Meters::new(2400.0),
};
