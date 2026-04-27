// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory Catalog Module
//!
//! This module provides constant definitions for well-known observatories.
//!
//! Named observatories are represented as [`Observatory`] values, carrying a
//! geodetic location and a reference atmospheric pressure.  Unnamed or
//! pressure-unknown sites are still available as bare
//! [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic) constants.
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
//! ```rust
//! use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
//! use siderust::qtty::Meter;
//!
//! let ecef = ROQUE_DE_LOS_MUCHACHOS.geodetic().to_cartesian::<Meter>();
//! let p_hpa = ROQUE_DE_LOS_MUCHACHOS.reference_pressure.value();
//! ```

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::{Degrees, Hectopascals, Meters};

/// A well-known astronomical observatory with location and reference atmosphere.
pub struct Observatory {
    /// Human-readable name.
    pub name: &'static str,
    /// Geodetic position (WGS84 ellipsoidal lon/lat/height).
    pub geodetic: Geodetic<ECEF>,
    /// Median reference atmospheric pressure for the site.
    pub reference_pressure: Hectopascals,
}

impl Observatory {
    /// Returns the geodetic position of this observatory.
    #[inline]
    pub const fn geodetic(&self) -> Geodetic<ECEF> {
        self.geodetic
    }
}

/// Roque de los Muchachos Observatory (La Palma, Canary Islands, Spain).
///
/// - Longitude: −17.8925°
/// - Latitude:  +28.7543°
/// - Altitude:  2396 m (WGS84 ellipsoidal height)
/// - Reference pressure: 744 hPa (source: `NSB_Utils.py:59`, `la_palma_pres = 744 hPa`)
///
/// For an atmosphere profile calibrated to this site, see
/// [`crate::atmosphere::profile::AtmosphereProfile::LA_PALMA`].
pub const ROQUE_DE_LOS_MUCHACHOS: Observatory = Observatory {
    name: "Roque de los Muchachos Observatory",
    geodetic: Geodetic::new_raw(
        Degrees::new(-17.8925),
        Degrees::new(28.7543),
        Meters::new(2396.0),
    ),
    reference_pressure: Hectopascals::new(744.0),
};

/// El Paranal Observatory (Atacama Desert, Chile).
///
/// - Longitude: −70.4043°
/// - Latitude:  −24.6272°
/// - Altitude:  2635 m (WGS84 ellipsoidal height)
/// - Reference pressure: 744 hPa (source: `NSB_Utils.py:59`, `cerro_paranal_pres = 744 hPa`)
pub const EL_PARANAL: Observatory = Observatory {
    name: "El Paranal Observatory",
    geodetic: Geodetic::new_raw(
        Degrees::new(-70.4043),
        Degrees::new(-24.6272),
        Meters::new(2635.0),
    ),
    reference_pressure: Hectopascals::new(744.0),
};

/// Mauna Kea Observatory (Hawaiʻi, USA)
///
/// - Longitude: −155.4681°
/// - Latitude:  +19.8207°
/// - Altitude:  4207 m (WGS84 ellipsoidal height)
///
/// For an atmosphere profile calibrated to this site, see
/// [`crate::atmosphere::profile::AtmosphereProfile::MAUNA_KEA`].
pub const MAUNA_KEA: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-155.4681),
    Degrees::new(19.8207),
    Meters::new(4207.0),
);

/// La Silla Observatory (ESO, Chile)
///
/// - Longitude: −70.7346°
/// - Latitude:  −29.2584°
/// - Altitude:  2400 m (WGS84 ellipsoidal height)
///
/// For an atmosphere profile calibrated to this site, see
/// [`crate::atmosphere::profile::AtmosphereProfile::LA_SILLA`].
pub const LA_SILLA_OBSERVATORY: Geodetic<ECEF> = Geodetic::new_raw(
    Degrees::new(-70.7346),
    Degrees::new(-29.2584),
    Meters::new(2400.0),
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn el_paranal_reference_pressure() {
        assert_eq!(EL_PARANAL.reference_pressure.value(), 744.0,
            "El Paranal reference pressure should be 744 hPa (NSB_Utils.py:59)");
    }

    #[test]
    fn roque_de_los_muchachos_reference_pressure() {
        assert_eq!(ROQUE_DE_LOS_MUCHACHOS.reference_pressure.value(), 744.0,
            "Roque de los Muchachos reference pressure should be 744 hPa (NSB_Utils.py:59)");
    }
}
