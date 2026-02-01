//! Frame-specific Position constructors with astronomical conventions.

use super::{clamp_polar, normalize_azimuth};
use super::position_core::Position;
use crate::coordinates::{centers, frames};
use qtty::{Degrees, LengthUnit, Quantity};

// =============================================================================
// ICRS
// =============================================================================

impl<C, U> Position<C, frames::ICRS, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates an ICRS position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
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

impl<C, U> Position<C, frames::EquatorialMeanJ2000, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a mean-J2000 equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
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

impl<C, U> Position<C, frames::EquatorialMeanOfDate, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a mean-of-date equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
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

impl<C, U> Position<C, frames::EquatorialTrueOfDate, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a true-of-date equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
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

impl<C, U> Position<C, frames::Ecliptic, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates an Ecliptic position from longitude, latitude, and distance.
    ///
    /// # Arguments
    /// - `lon`: Ecliptic longitude (λ), normalized to [0°, 360°)
    /// - `lat`: Ecliptic latitude (β), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(lon: Degrees, lat: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(lat),
            normalize_azimuth(lon),
            distance.into(),
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

impl<C, U> Position<C, frames::Horizontal, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a Horizontal position from altitude, azimuth, and distance (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `alt`: Altitude (elevation), clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    /// - `distance`: Radial distance from the observer
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(alt: Degrees, az: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(alt),
            normalize_azimuth(az),
            distance.into(),
        ))
    }

    /// Returns the altitude in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

impl<U: LengthUnit> Position<centers::Topocentric, frames::Horizontal, U> {
    /// Creates a Topocentric Horizontal position with observer site parameters (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `site`: The observer site parameters
    /// - `alt`: Altitude (elevation), clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    /// - `distance`: Radial distance from the observer
    #[inline]
    pub fn new_with_site<T: Into<Quantity<U>>>(
        site: centers::ObserverSite,
        alt: Degrees,
        az: Degrees,
        distance: T,
    ) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw_with_params(
            site,
            clamp_polar(alt),
            normalize_azimuth(az),
            distance.into(),
        ))
    }

    /// Returns the altitude in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

// =============================================================================
// Geographic (ECEF)
// =============================================================================

impl<C, U> Position<C, frames::ECEF, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a Geographic position from longitude, latitude, and altitude.
    ///
    /// # Arguments
    /// - `lon`: Geodetic longitude, positive eastward, normalized to [0°, 360°)
    /// - `lat`: Geodetic latitude, positive northward, clamped to [-90°, +90°]
    /// - `altitude`: Height above the WGS84 ellipsoid (e.g., observatory altitude)
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(lon: Degrees, lat: Degrees, altitude: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(lat),
            normalize_azimuth(lon),
            altitude.into(),
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
