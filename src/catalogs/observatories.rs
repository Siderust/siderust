// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory catalog
//!
//! ## Scientific scope
//!
//! Ground-based astronomical observations are referenced to a *site*: a
//! geodetic position on the WGS84 ellipsoid plus a typical local
//! atmospheric state (pressure, temperature, relative humidity). These
//! quantities feed every site-dependent computation in the crate — airmass,
//! refraction, atmospheric extinction, night-sky brightness, parallactic
//! correction — so canonical, citation-backed constants for the major
//! observatories live here in one place.
//!
//! Positions are stored as ellipsoidal **(longitude, latitude, ellipsoidal
//! height above WGS84)** triples, *not* as spherical positions. A
//! geodetic latitude is the angle between the local vertical and the
//! equatorial plane on the WGS84 reference ellipsoid; it does not coincide
//! with a geocentric latitude, and the ellipsoidal height is **not** a
//! radial distance. Converting to ECEF therefore requires the
//! ellipsoid-aware Bowring formula (see
//! [`affn::ellipsoidal::Position::to_cartesian`]) rather than a naïve
//! spherical-to-Cartesian unrolling.
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - [`Observatory`] — site record carrying name, [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic)
//!   position, [`Hectopascals`] reference pressure, optional [`Kelvins`]
//!   reference temperature, and an optional dimensionless reference
//!   relative humidity in `[0.0, 1.0]`.
//! - Named constants: [`ROQUE_DE_LOS_MUCHACHOS`], [`EL_PARANAL`],
//!   [`MAUNA_KEA`], [`LA_SILLA_OBSERVATORY`].
//!
//! Reference relative humidity is intentionally a bare `f64` (no
//! `RelativeHumidity` newtype exists yet); it is documented as a
//! dimensionless fraction in `[0.0, 1.0]` (so 14.5 % RH is `0.145`, never
//! `14.5`).
//!
//! Atmospheric profiles built from these site constants live in the
//! `atmosphere` module (behind the `atmosphere` feature).
//!
//! ## Usage
//!
//! ```rust
//! use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
//! use siderust::qtty::Meter;
//!
//! let ecef = ROQUE_DE_LOS_MUCHACHOS.geodetic().to_cartesian::<Meter>();
//! let p_hpa = ROQUE_DE_LOS_MUCHACHOS.reference_pressure.value();
//! ```
//!
//! ## References
//!
//! - National Imagery and Mapping Agency (2000). *Department of Defense
//!   World Geodetic System 1984: Its Definition and Relationships with
//!   Local Geodetic Systems*. NIMA Technical Report TR8350.2, 3rd ed.
//! - Bowring, B. R. (1976). "Transformation from spatial to geographical
//!   coordinates". *Survey Review* **23**(181), 323–327.
//! - Patat, F., Moehler, S., O'Brien, K., et al. (2011). "Optical
//!   atmospheric extinction over Cerro Paranal". *Astronomy &
//!   Astrophysics* **527**, A91. doi:10.1051/0004-6361/201015537.
//! - Instituto de Astrofísica de Canarias. *Site characterization of the
//!   Observatorio del Roque de los Muchachos*. <https://www.iac.es/en/observatorios-de-canarias>.

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
    /// (e.g. `AtmosphereProfile` from the `atmosphere` feature) rather than
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
/// - Reference relative humidity: `None` — same rationale as temperature.
///
/// For an atmosphere profile calibrated to this site, see
/// `AtmosphereProfile::ROQUE_DE_LOS_MUCHACHOS` (requires the `atmosphere` feature).
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

/// Mauna Kea Observatory (Hawaiʻi, USA).
///
/// - Longitude: −155.4681°
/// - Latitude:  +19.8207°
/// - Altitude:  4207 m (WGS84 ellipsoidal height)
/// - Reference pressure: 614 hPa (ISA barometric formula at 4207 m altitude)
/// - Reference temperature: `None` — no authoritative cited value wired in yet.
/// - Reference relative humidity: `None` — same rationale as temperature.
///
/// For an atmosphere profile calibrated to this site, see
/// `AtmosphereProfile::MAUNA_KEA` (requires the `atmosphere` feature).
pub const MAUNA_KEA: Observatory = Observatory {
    name: "Mauna Kea Observatory",
    geodetic: Geodetic::new_raw(
        Degrees::new(-155.4681),
        Degrees::new(19.8207),
        Meters::new(4207.0),
    ),
    reference_pressure: Hectopascals::new(614.0),
    reference_temperature: None,
    reference_relative_humidity: None,
};

/// La Silla Observatory (ESO, Chile).
///
/// - Longitude: −70.7346°
/// - Latitude:  −29.2584°
/// - Altitude:  2400 m (WGS84 ellipsoidal height)
/// - Reference pressure: 744 hPa (same altitude class as El Paranal and Roque)
/// - Reference temperature: `None` — no authoritative cited value wired in yet.
/// - Reference relative humidity: `None` — same rationale as temperature.
///
/// For an atmosphere profile calibrated to this site, see
/// `AtmosphereProfile::LA_SILLA` (requires the `atmosphere` feature).
pub const LA_SILLA_OBSERVATORY: Observatory = Observatory {
    name: "La Silla Observatory",
    geodetic: Geodetic::new_raw(
        Degrees::new(-70.7346),
        Degrees::new(-29.2584),
        Meters::new(2400.0),
    ),
    reference_pressure: Hectopascals::new(744.0),
    reference_temperature: None,
    reference_relative_humidity: None,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn el_paranal_reference_pressure() {
        assert_eq!(
            EL_PARANAL.reference_pressure.value(),
            744.0,
            "El Paranal reference pressure should be 744 hPa (NSB_Utils.py:59)"
        );
    }

    #[test]
    fn roque_de_los_muchachos_reference_pressure() {
        assert_eq!(
            ROQUE_DE_LOS_MUCHACHOS.reference_pressure.value(),
            744.0,
            "Roque de los Muchachos reference pressure should be 744 hPa (NSB_Utils.py:59)"
        );
    }

    #[test]
    fn mauna_kea_reference_pressure() {
        assert_eq!(
            MAUNA_KEA.reference_pressure.value(),
            614.0,
            "Mauna Kea reference pressure should be ~614 hPa at 4207 m altitude"
        );
    }

    #[test]
    fn la_silla_reference_pressure() {
        assert_eq!(
            LA_SILLA_OBSERVATORY.reference_pressure.value(),
            744.0,
            "La Silla reference pressure should be 744 hPa"
        );
    }

    #[test]
    fn el_paranal_optional_atmosphere_unset() {
        assert!(
            EL_PARANAL.reference_temperature.is_none(),
            "El Paranal reference temperature must stay None until a cited value is wired in"
        );
        assert!(
            EL_PARANAL.reference_relative_humidity.is_none(),
            "El Paranal reference relative humidity must stay None until a cited value is wired in"
        );
    }

    #[test]
    fn roque_de_los_muchachos_optional_atmosphere_unset() {
        assert!(
            ROQUE_DE_LOS_MUCHACHOS.reference_temperature.is_none(),
            "ORM reference temperature must stay None until a cited value is wired in"
        );
        assert!(
            ROQUE_DE_LOS_MUCHACHOS.reference_relative_humidity.is_none(),
            "ORM reference relative humidity must stay None until a cited value is wired in"
        );
    }

    #[test]
    fn mauna_kea_optional_atmosphere_unset() {
        assert!(MAUNA_KEA.reference_temperature.is_none());
        assert!(MAUNA_KEA.reference_relative_humidity.is_none());
    }

    #[test]
    fn la_silla_optional_atmosphere_unset() {
        assert!(LA_SILLA_OBSERVATORY.reference_temperature.is_none());
        assert!(LA_SILLA_OBSERVATORY.reference_relative_humidity.is_none());
    }
}
