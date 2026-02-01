//! Frame-specific Direction constructors with astronomical conventions.

use super::{clamp_polar, normalize_azimuth};
use super::direction_core::Direction;
use crate::coordinates::frames;
use qtty::Degrees;

// =============================================================================
// ICRS
// =============================================================================

impl Direction<frames::ICRS> {
    /// Creates an ICRS direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Equatorial Frames
// =============================================================================

impl Direction<frames::EquatorialMeanJ2000> {
    /// Creates a mean-J2000 equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl Direction<frames::EquatorialMeanOfDate> {
    /// Creates a mean-of-date equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl Direction<frames::EquatorialTrueOfDate> {
    /// Creates a true-of-date equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Ecliptic
// =============================================================================

impl Direction<frames::Ecliptic> {
    /// Creates an Ecliptic direction from ecliptic longitude and latitude.
    ///
    /// # Arguments
    /// - `lon`: Ecliptic longitude (λ), normalized to [0°, 360°)
    /// - `lat`: Ecliptic latitude (β), clamped to [-90°, +90°]
    #[inline]
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(lat),
            normalize_azimuth(lon),
        ))
    }

    /// Returns the ecliptic longitude (λ) in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the ecliptic latitude (β) in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Horizontal
// =============================================================================

impl Direction<frames::Horizontal> {
    /// Creates a Horizontal direction from altitude and azimuth (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `alt`: Altitude (elevation) above the horizon, clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    #[inline]
    pub fn new(alt: Degrees, az: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(alt),
            normalize_azimuth(az),
        ))
    }

    /// Returns the altitude (elevation) in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth (from North through East) in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

// =============================================================================
// Geographic (ECEF)
// =============================================================================

impl Direction<frames::ECEF> {
    /// Creates a Geographic direction from longitude and latitude.
    ///
    /// # Arguments
    /// - `lon`: Geodetic longitude, positive eastward, normalized to [0°, 360°)
    /// - `lat`: Geodetic latitude, positive northward, clamped to [-90°, +90°]
    #[inline]
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(lat),
            normalize_azimuth(lon),
        ))
    }

    /// Returns the geodetic longitude in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the geodetic latitude in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}
