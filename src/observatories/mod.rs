//! # Observatory Catalog Module
//!
//! This module provides a set of constant definitions for well-known observatories,
//! represented by their geographic position. (see [`crate::coordinates::spherical::GeographicCoord`]).

use crate::coordinates::spherical::GeographicCoord;
use crate::units::{Degrees, Kilometers};

/// Roque de los Muchachos Observatory (La Palma, Canary Islands, Spain)
/// - Longitude: −17.8925°
/// - Latitude:  +28.7543°
/// - Altitude:  2396 m
pub const ROQUE_DE_LOS_MUCHACHOS: GeographicCoord = GeographicCoord::new_const(
    Degrees::new(-17.8925),
    Degrees::new(28.7543),
    Some(Kilometers::new(2.396)),
);

/// El Panaral (Chile)
/// - Longitude: −70.4043°
/// - Latitude:  −24.6272°
/// - Altitude:  2635 m
pub const EL_PANARAL: GeographicCoord = GeographicCoord::new_const(
    Degrees::new(-70.4043),
    Degrees::new(-24.6272),
    Some(Kilometers::new(2.635)),
);

/// Mauna Kea Observatory (Hawaiʻi, USA)
/// - Longitude: −155.4681°
/// - Latitude:  +19.8207°
/// - Altitude:  4207 m
pub const MAUNA_KEA: GeographicCoord = GeographicCoord::new_const(
    Degrees::new(-155.4681),
    Degrees::new(19.8207),
    Some(Kilometers::new(4.207)),
);

/// La Silla Observatory (ESO, Chile)
/// - Longitude: −70.7346°
/// - Latitude:  −29.2584°
/// - Altitude:  2400 m
pub const LA_SILLA_OBSERVATORY: GeographicCoord = GeographicCoord::new_const(
    Degrees::new(-70.7346),
    Degrees::new(-29.2584),
    Some(Kilometers::new(2.400)),
);
