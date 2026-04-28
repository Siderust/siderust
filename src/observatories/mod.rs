// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory Catalog Module
//!
//! This module provides constant definitions for well-known observatories.
//!
//! Named observatories are represented as [`Observatory`] values, carrying a
//! geodetic location, a reference atmospheric pressure, and optional reference
//! temperature and relative-humidity values.  Unnamed or pressure-unknown sites
//! are still available as bare
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
use crate::qtty::{Degrees, Hectopascals, Kelvins, Meters};

/// A well-known astronomical observatory with location and reference atmosphere.
pub struct Observatory {
    /// Human-readable name.
    pub name: &'static str,
    /// Geodetic position (WGS84 ellipsoidal lon/lat/height).
    pub geodetic: Geodetic<ECEF>,
    /// Median reference atmospheric pressure for the site.
    pub reference_pressure: Hectopascals,
    /// Optional median reference air temperature for the site.
    ///
    /// Always expressed on the absolute thermodynamic (kelvin) scale; convert
    /// any Celsius/Fahrenheit observation to kelvin at the boundary before
    /// constructing this value.  `None` means "no authoritative cited value
    /// has been wired in yet" — callers should fall back to a profile default
    /// (e.g. [`crate::atmosphere::profile::AtmosphereProfile`]) rather than
    /// assume a hard-coded constant.
    pub reference_temperature: Option<Kelvins>,
    /// Optional median reference relative humidity for the site.
    ///
    /// Encoded as a dimensionless fraction in `[0.0, 1.0]` (so 14.5 % RH is
    /// `0.145`, not `14.5`).  `None` carries the same "no authoritative cited
    /// value yet" semantics as [`Observatory::reference_temperature`].
    pub reference_relative_humidity: Option<f64>,
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
/// - Reference temperature: `None` — no authoritative cited value wired in yet.
///   The upstream NSB code does not carry a per-site temperature constant; an
///   ORM site-characterization citation (IAC) should be added in a follow-up.
/// - Reference relative humidity: `None` — same rationale as temperature.
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
    reference_temperature: None,
    reference_relative_humidity: None,
};

/// El Paranal Observatory (Atacama Desert, Chile).
///
/// - Longitude: −70.4043°
/// - Latitude:  −24.6272°
/// - Altitude:  2635 m (WGS84 ellipsoidal height)
/// - Reference pressure: 744 hPa (source: `NSB_Utils.py:59`, `cerro_paranal_pres = 744 hPa`)
/// - Reference temperature: `None` — no authoritative cited value wired in yet.
///   The upstream NSB code does not carry a per-site temperature constant; an
///   ESO Paranal site-characterization citation (e.g. Patat et al. 2011) should
///   be added in a follow-up.
/// - Reference relative humidity: `None` — same rationale as temperature.
pub const EL_PARANAL: Observatory = Observatory {
    name: "El Paranal Observatory",
    geodetic: Geodetic::new_raw(
        Degrees::new(-70.4043),
        Degrees::new(-24.6272),
        Meters::new(2635.0),
    ),
    reference_pressure: Hectopascals::new(744.0),
    reference_temperature: None,
    reference_relative_humidity: None,
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

    #[test]
    fn el_paranal_optional_atmosphere_unset() {
        assert!(EL_PARANAL.reference_temperature.is_none(),
            "El Paranal reference temperature must stay None until a cited value is wired in");
        assert!(EL_PARANAL.reference_relative_humidity.is_none(),
            "El Paranal reference relative humidity must stay None until a cited value is wired in");
    }

    #[test]
    fn roque_de_los_muchachos_optional_atmosphere_unset() {
        assert!(ROQUE_DE_LOS_MUCHACHOS.reference_temperature.is_none(),
            "ORM reference temperature must stay None until a cited value is wired in");
        assert!(ROQUE_DE_LOS_MUCHACHOS.reference_relative_humidity.is_none(),
            "ORM reference relative humidity must stay None until a cited value is wired in");
    }
}
