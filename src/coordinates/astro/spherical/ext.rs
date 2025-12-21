//! Extension traits for frame-specific spherical coordinate methods.
//!
//! These traits provide astronomy-specific convenience methods on `affn` coordinate
//! types without violating Rust's orphan rules.
//!
//! # Design
//!
//! Instead of adding inherent methods like:
//! ```ignore
//! impl Direction<ICRS> {
//!     pub fn new_icrs(ra: Degrees, dec: Degrees) -> Self { ... }
//! }
//! ```
//!
//! We use extension traits:
//! ```ignore
//! pub trait IcrsDirectionExt {
//!     fn new_icrs(ra: Degrees, dec: Degrees) -> Self;
//!     fn ra(&self) -> Degrees;
//!     fn dec(&self) -> Degrees;
//! }
//! ```
//!
//! This is zero-cost (all methods are inlined) and respects coherence rules.

use crate::coordinates::algebra::centers::ReferenceCenter;
use crate::coordinates::algebra::frames::{Ecliptic, Equatorial, Horizontal, ECEF, ICRS};
use crate::coordinates::algebra::spherical::{Direction, Position};
use qtty::{Degrees, LengthUnit, Quantity};

// =============================================================================
// ICRS Extension Traits
// =============================================================================

/// Extension trait for ICRS spherical directions.
///
/// Provides astronomy-specific constructors and accessors using
/// Right Ascension (RA) and Declination (Dec) terminology.
pub trait IcrsDirectionExt {
    /// Creates a new ICRS direction with RA and Dec.
    ///
    /// Right Ascension is normalized to [0°, 360°], Declination to [-90°, 90°].
    fn new_icrs(ra: Degrees, dec: Degrees) -> Self;
    
    /// Returns the Right Ascension (α) in degrees.
    fn ra(&self) -> Degrees;
    
    /// Returns the Declination (δ) in degrees.
    fn dec(&self) -> Degrees;
}

impl IcrsDirectionExt for Direction<ICRS> {
    #[inline]
    fn new_icrs(ra: Degrees, dec: Degrees) -> Self {
        Self::new(dec.wrap_quarter_fold(), ra.normalize())
    }
    
    #[inline]
    fn ra(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn dec(&self) -> Degrees {
        self.polar
    }
}

/// Extension trait for ICRS spherical positions.
pub trait IcrsPositionExt<U: LengthUnit> {
    /// Creates a new ICRS position with RA, Dec, and distance.
    fn new_icrs<A: Into<Degrees>, T: Into<Quantity<U>>>(ra: A, dec: A, distance: T) -> Self;
    
    /// Returns the Right Ascension (α) in degrees.
    fn ra(&self) -> Degrees;
    
    /// Returns the Declination (δ) in degrees.
    fn dec(&self) -> Degrees;
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> IcrsPositionExt<U> for Position<C, ICRS, U> {
    #[inline]
    fn new_icrs<A: Into<Degrees>, T: Into<Quantity<U>>>(ra: A, dec: A, distance: T) -> Self {
        Self::new_raw(
            dec.into().wrap_quarter_fold(),
            ra.into().normalize(),
            distance.into(),
        )
    }
    
    #[inline]
    fn ra(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn dec(&self) -> Degrees {
        self.polar
    }
}

// =============================================================================
// Equatorial Extension Traits
// =============================================================================

/// Extension trait for Equatorial spherical directions.
pub trait EquatorialDirectionExt {
    /// Creates a new Equatorial direction with RA and Dec.
    fn new_equatorial(ra: Degrees, dec: Degrees) -> Self;
    
    /// Returns the Right Ascension (α) in degrees.
    fn ra(&self) -> Degrees;
    
    /// Returns the Declination (δ) in degrees.
    fn dec(&self) -> Degrees;
}

impl EquatorialDirectionExt for Direction<Equatorial> {
    #[inline]
    fn new_equatorial(ra: Degrees, dec: Degrees) -> Self {
        Self::new(dec.wrap_quarter_fold(), ra.normalize())
    }
    
    #[inline]
    fn ra(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn dec(&self) -> Degrees {
        self.polar
    }
}

/// Extension trait for Equatorial spherical positions.
pub trait EquatorialPositionExt<U: LengthUnit> {
    /// Creates a new Equatorial position with RA, Dec, and distance.
    fn new_equatorial<A: Into<Degrees>, T: Into<Quantity<U>>>(ra: A, dec: A, distance: T) -> Self;
    
    /// Returns the Right Ascension (α) in degrees.
    fn ra(&self) -> Degrees;
    
    /// Returns the Declination (δ) in degrees.
    fn dec(&self) -> Degrees;
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> EquatorialPositionExt<U>
    for Position<C, Equatorial, U>
{
    #[inline]
    fn new_equatorial<A: Into<Degrees>, T: Into<Quantity<U>>>(ra: A, dec: A, distance: T) -> Self {
        Self::new_raw(
            dec.into().wrap_quarter_fold(),
            ra.into().normalize(),
            distance.into(),
        )
    }
    
    #[inline]
    fn ra(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn dec(&self) -> Degrees {
        self.polar
    }
}

// =============================================================================
// Ecliptic Extension Traits
// =============================================================================

/// Extension trait for Ecliptic spherical directions.
pub trait EclipticDirectionExt {
    /// Creates a new Ecliptic direction with longitude and latitude.
    fn new_ecliptic(lon: Degrees, lat: Degrees) -> Self;
    
    /// Returns the ecliptic longitude (λ) in degrees.
    fn lon(&self) -> Degrees;
    
    /// Returns the ecliptic latitude (β) in degrees.
    fn lat(&self) -> Degrees;
}

impl EclipticDirectionExt for Direction<Ecliptic> {
    #[inline]
    fn new_ecliptic(lon: Degrees, lat: Degrees) -> Self {
        Self::new(lat.wrap_quarter_fold(), lon.normalize())
    }
    
    #[inline]
    fn lon(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn lat(&self) -> Degrees {
        self.polar
    }
}

/// Extension trait for Ecliptic spherical positions.
pub trait EclipticPositionExt<U: LengthUnit> {
    /// Creates a new Ecliptic position with longitude, latitude, and distance.
    fn new_ecliptic<A: Into<Degrees>, T: Into<Quantity<U>>>(lon: A, lat: A, distance: T) -> Self;
    
    /// Returns the ecliptic longitude (λ) in degrees.
    fn lon(&self) -> Degrees;
    
    /// Returns the ecliptic latitude (β) in degrees.
    fn lat(&self) -> Degrees;
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> EclipticPositionExt<U>
    for Position<C, Ecliptic, U>
{
    #[inline]
    fn new_ecliptic<A: Into<Degrees>, T: Into<Quantity<U>>>(lon: A, lat: A, distance: T) -> Self {
        Self::new_raw(
            lat.into().wrap_quarter_fold(),
            lon.into().normalize(),
            distance.into(),
        )
    }
    
    #[inline]
    fn lon(&self) -> Degrees {
        self.azimuth
    }
    
    #[inline]
    fn lat(&self) -> Degrees {
        self.polar
    }
}

// =============================================================================
// Horizontal Extension Traits
// =============================================================================

/// Extension trait for Horizontal spherical directions.
pub trait HorizontalDirectionExt {
    /// Creates a new Horizontal direction with altitude and azimuth.
    fn new_horizontal(alt: Degrees, az: Degrees) -> Self;
    
    /// Returns the altitude (elevation) in degrees.
    fn alt(&self) -> Degrees;
    
    /// Returns the azimuth in degrees.
    fn az(&self) -> Degrees;
}

impl HorizontalDirectionExt for Direction<Horizontal> {
    #[inline]
    fn new_horizontal(alt: Degrees, az: Degrees) -> Self {
        Self::new(alt.wrap_quarter_fold(), az.normalize())
    }
    
    #[inline]
    fn alt(&self) -> Degrees {
        self.polar
    }
    
    #[inline]
    fn az(&self) -> Degrees {
        self.azimuth
    }
}

/// Extension trait for Horizontal spherical positions.
pub trait HorizontalPositionExt<U: LengthUnit> {
    /// Creates a new Horizontal position with altitude, azimuth, and distance.
    /// Only available for centers with no parameters (e.g., Geocentric).
    fn new_horizontal<A: Into<Degrees>, T: Into<Quantity<U>>>(alt: A, az: A, distance: T) -> Self;
    
    /// Returns the altitude (elevation) in degrees.
    fn alt(&self) -> Degrees;
    
    /// Returns the azimuth in degrees.
    fn az(&self) -> Degrees;
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> HorizontalPositionExt<U>
    for Position<C, Horizontal, U>
{
    #[inline]
    fn new_horizontal<A: Into<Degrees>, T: Into<Quantity<U>>>(alt: A, az: A, distance: T) -> Self {
        Self::new_raw(
            alt.into().wrap_quarter_fold(),
            az.into().normalize(),
            distance.into(),
        )
    }
    
    #[inline]
    fn alt(&self) -> Degrees {
        self.polar
    }
    
    #[inline]
    fn az(&self) -> Degrees {
        self.azimuth
    }
}

/// Extension trait for reading Horizontal position coordinates (works with Topocentric center).
pub trait HorizontalPositionReadExt {
    /// Returns the altitude (elevation) in degrees.
    fn alt(&self) -> Degrees;
    
    /// Returns the azimuth in degrees.
    fn az(&self) -> Degrees;
}

impl<U: LengthUnit> HorizontalPositionReadExt
    for Position<crate::coordinates::centers::Topocentric, Horizontal, U>
{
    #[inline]
    fn alt(&self) -> Degrees {
        self.polar
    }
    
    #[inline]
    fn az(&self) -> Degrees {
        self.azimuth
    }
}

// =============================================================================
// ECEF (Geographic) Extension Traits
// =============================================================================

/// Extension trait for ECEF/Geographic spherical directions.
pub trait EcefDirectionExt {
    /// Creates a new geographic direction with latitude and longitude.
    fn new_geographic(lat: Degrees, lon: Degrees) -> Self;
    
    /// Returns the geodetic latitude (φ) in degrees.
    fn lat(&self) -> Degrees;
    
    /// Returns the geodetic longitude (λ) in degrees.
    fn lon(&self) -> Degrees;
}

impl EcefDirectionExt for Direction<ECEF> {
    #[inline]
    fn new_geographic(lat: Degrees, lon: Degrees) -> Self {
        Self::new(lat.wrap_quarter_fold(), lon.normalize())
    }
    
    #[inline]
    fn lat(&self) -> Degrees {
        self.polar
    }
    
    #[inline]
    fn lon(&self) -> Degrees {
        self.azimuth
    }
}

/// Extension trait for ECEF/Geographic spherical positions.
pub trait EcefPositionExt<U: LengthUnit> {
    /// Creates a new geographic position with latitude, longitude, and distance.
    fn new_geographic<A: Into<Degrees>, T: Into<Quantity<U>>>(lat: A, lon: A, distance: T) -> Self;
    
    /// Returns the geodetic latitude (φ) in degrees.
    fn lat(&self) -> Degrees;
    
    /// Returns the geodetic longitude (λ) in degrees.
    fn lon(&self) -> Degrees;
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> EcefPositionExt<U> for Position<C, ECEF, U> {
    #[inline]
    fn new_geographic<A: Into<Degrees>, T: Into<Quantity<U>>>(lat: A, lon: A, distance: T) -> Self {
        Self::new_raw(
            lat.into().wrap_quarter_fold(),
            lon.into().normalize(),
            distance.into(),
        )
    }
    
    #[inline]
    fn lat(&self) -> Degrees {
        self.polar
    }
    
    #[inline]
    fn lon(&self) -> Degrees {
        self.azimuth
    }
}
