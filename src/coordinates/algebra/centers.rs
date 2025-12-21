//! # Reference Centers Module
//!
//! This module defines astronomical reference centers (origins) for coordinate systems.
//! A reference center specifies the origin point from which positions are measured.
//!
//! ## Architecture
//!
//! All center types implement the [`ReferenceCenter`] trait from `affn`, which provides
//! a common interface. The trait itself is re-exported from `affn` for convenience.
//!
//! ## Predefined Centers
//!
//! The following reference centers are provided:
//!
//! - [`Barycentric`]: Center of mass of the solar system.
//! - [`Heliocentric`]: Center of the Sun.
//! - [`Geocentric`]: Center of the Earth.
//! - [`Topocentric`]: Observer's location on the surface of the Earth (parameterized by [`ObserverSite`]).
//! - [`Bodycentric`]: Generic center for any orbiting celestial body (parameterized by [`BodycentricParams`]).
//!
//! ## Parameterized Centers
//!
//! Some reference centers require runtime parameters:
//!
//! - For most centers (Barycentric, Heliocentric, Geocentric), `Params = ()` (zero-cost).
//! - For [`Topocentric`], `Params = ObserverSite` which stores the observer's geographic location.
//! - For [`Bodycentric`], `Params = BodycentricParams` which stores the body's orbital elements.
//!
//! ## Extending
//!
//! To define a new reference center, use the [`affn::new_center!`] macro:
//!
//! ```rust
//! use affn::new_center;
//! use affn::ReferenceCenter;
//!
//! new_center!(Lunarcentric);
//! assert_eq!(Lunarcentric::center_name(), "Lunarcentric");
//! ```
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::centers::{ReferenceCenter, Geocentric};
//!
//! let name = Geocentric::center_name();
//! assert_eq!(name, "Geocentric");
//! ```

use crate::astro::orbit::Orbit;
use qtty::{Degrees, Meter, Quantity};
use std::fmt::Debug;

// Re-export core traits from affn
pub use affn::{AffineCenter, NoCenter, ReferenceCenter};

// Required for Transform specialization
#[derive(Debug, Copy, Clone)]
pub struct Heliocentric;

impl ReferenceCenter for Heliocentric {
    type Params = ();

    fn center_name() -> &'static str {
        stringify!(Heliocentric)
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Barycentric;

impl ReferenceCenter for Barycentric {
    type Params = ();

    fn center_name() -> &'static str {
        stringify!(Barycentric)
    }
}

// =============================================================================
// ObserverSite: Parameters for Topocentric coordinates
// =============================================================================

/// Geographic location of an observer, used as parameters for [`Topocentric`] coordinates.
///
/// This struct stores the observer's position on or above Earth's surface.
/// It is embedded in coordinate values when the center is [`Topocentric`],
/// allowing transformations to use the site information without external context.
///
/// # Fields
///
/// - `lon`: Geodetic longitude, positive eastward, in degrees.
/// - `lat`: Geodetic latitude, positive northward, in degrees.
/// - `height`: Height above the reference ellipsoid (WGS84), in meters.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::centers::ObserverSite;
/// use qtty::*;
///
/// let greenwich = ObserverSite {
///     lon: 0.0 * DEG,
///     lat: 51.4769 * DEG,
///     height: 0.0 * M,
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverSite {
    /// Geodetic longitude (positive eastward), in degrees.
    pub lon: Degrees,
    /// Geodetic latitude (positive northward), in degrees.
    pub lat: Degrees,
    /// Height above the WGS84 ellipsoid, in meters.
    pub height: Quantity<Meter>,
}

impl Default for ObserverSite {
    /// Returns an observer at the origin (0°, 0°, 0m).
    ///
    /// Note: This default is primarily for internal use. In practice, you should
    /// always provide meaningful site coordinates for topocentric calculations.
    fn default() -> Self {
        Self {
            lon: Degrees::new(0.0),
            lat: Degrees::new(0.0),
            height: Quantity::<Meter>::new(0.0),
        }
    }
}

impl ObserverSite {
    /// Creates a new observer site from longitude, latitude, and height.
    ///
    /// # Arguments
    ///
    /// - `lon`: Geodetic longitude (positive eastward), in degrees.
    /// - `lat`: Geodetic latitude (positive northward), in degrees.
    /// - `height`: Height above the WGS84 ellipsoid, in meters.
    pub fn new(lon: Degrees, lat: Degrees, height: Quantity<Meter>) -> Self {
        Self { lon, lat, height }
    }

    /// Creates an `ObserverSite` from a [`Geographic`] spherical position.
    ///
    /// This is a convenience method for converting observatory locations
    /// (typically defined as `Geographic` positions) into the parameterized
    /// site representation used by `Topocentric` coordinates.
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::algebra::centers::ObserverSite;
    /// use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
    ///
    /// let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    /// ```
    pub fn from_geographic(geo: &crate::coordinates::astro::spherical::position::Geographic) -> Self {
        Self {
            lon: geo.azimuth, // longitude is stored in the azimuth field
            lat: geo.polar,   // latitude is stored in the polar field
            height: geo.distance.to::<Meter>(),
        }
    }

    /// Computes the observer's geocentric position in the ITRF/ECEF frame.
    ///
    /// This converts geodetic coordinates (latitude, longitude, height) to
    /// geocentric Cartesian coordinates using the WGS84 ellipsoid.
    ///
    /// # WGS84 Parameters
    ///
    /// - Semi-major axis (a): 6,378,137.0 m
    /// - Flattening (f): 1/298.257223563
    /// - First eccentricity squared (e²): 0.00669437999014
    ///
    /// # Returns
    ///
    /// A position vector in the ECEF (Earth-Centered Earth-Fixed) frame,
    /// with the specified length unit.
    ///
    /// # Note
    ///
    /// This position is in the ITRF/ECEF frame (rotating with Earth).
    /// To use it with celestial coordinates, you must rotate it to the
    /// appropriate celestial frame using Earth rotation parameters (GMST, etc.).
    pub fn geocentric_itrf<U: qtty::LengthUnit>(
        &self,
    ) -> crate::coordinates::algebra::cartesian::Position<Geocentric, crate::coordinates::algebra::frames::ECEF, U>
    where
        Quantity<U>: From<Quantity<Meter>>,
    {
        use qtty::Radian;

        // WGS84 ellipsoid parameters
        const A: f64 = 6_378_137.0; // Semi-major axis in meters
        const F: f64 = 1.0 / 298.257_223_563; // Flattening
        const E2: f64 = 2.0 * F - F * F; // First eccentricity squared

        // Convert geodetic to geocentric Cartesian (ECEF)
        let lat_rad = self.lat.to::<Radian>().value();
        let lon_rad = self.lon.to::<Radian>().value();
        let h = self.height.value();

        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_lon = lon_rad.sin();
        let cos_lon = lon_rad.cos();

        // Radius of curvature in the prime vertical
        let n = A / (1.0 - E2 * sin_lat * sin_lat).sqrt();

        // Geocentric Cartesian coordinates (meters)
        let x_m = (n + h) * cos_lat * cos_lon;
        let y_m = (n + h) * cos_lat * sin_lon;
        let z_m = (n * (1.0 - E2) + h) * sin_lat;

        // Convert to target units
        let x: Quantity<U> = Quantity::<Meter>::new(x_m).into();
        let y: Quantity<U> = Quantity::<Meter>::new(y_m).into();
        let z: Quantity<U> = Quantity::<Meter>::new(z_m).into();

        crate::coordinates::algebra::cartesian::Position::<Geocentric, crate::coordinates::algebra::frames::ECEF, U>::new(
            x, y, z
        )
    }
}

// =============================================================================
// Topocentric Center (parameterized)
// =============================================================================

/// Observer's location on the surface of the Earth.
///
/// Unlike other reference centers, `Topocentric` is *parameterized*: coordinates
/// with this center carry an [`ObserverSite`] that specifies the observer's
/// geographic location. This allows horizontal coordinates to know their
/// observation site without external context.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::algebra::centers::{Topocentric, ObserverSite, ReferenceCenter};
/// use qtty::*;
///
/// // Topocentric coordinates require an ObserverSite
/// let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
///
/// // The site is stored as Topocentric::Params
/// assert_eq!(std::mem::size_of::<<Topocentric as ReferenceCenter>::Params>(),
///            std::mem::size_of::<ObserverSite>());
/// ```
#[derive(Debug, Copy, Clone)]
pub struct Topocentric;

impl ReferenceCenter for Topocentric {
    type Params = ObserverSite;

    fn center_name() -> &'static str {
        "Topocentric"
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Geocentric;
impl ReferenceCenter for Geocentric {
    type Params = ();

    fn center_name() -> &'static str {
        "Geocentric"
    }
}

// =============================================================================
// Bodycentric: Generic center for any orbiting celestial body
// =============================================================================

/// Specifies the reference center for an orbit (where the orbit is defined relative to).
///
/// When transforming to/from body-centric coordinates, the orbit must be converted
/// to match the coordinate system being transformed. This enum indicates which
/// standard center the orbit is relative to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OrbitReferenceCenter {
    /// Orbit is defined relative to the solar system barycenter.
    Barycentric,
    /// Orbit is defined relative to the Sun (most common for planets, asteroids, comets).
    #[default]
    Heliocentric,
    /// Orbit is defined relative to Earth (for artificial satellites, the Moon).
    Geocentric,
}

/// Parameters for a body-centered coordinate system.
///
/// This struct specifies the orbital elements of a celestial body and which
/// reference center the orbit is defined relative to. This allows computing
/// the body's position at any Julian date using Keplerian propagation.
///
/// # Fields
///
/// - `orbit`: The Keplerian orbital elements of the body.
/// - `orbit_center`: Which center the orbit is defined relative to.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::algebra::centers::{BodycentricParams, OrbitReferenceCenter};
/// use siderust::astro::orbit::Orbit;
/// use siderust::astro::JulianDate;
/// use qtty::*;
///
/// // Mars-like orbit (heliocentric)
/// let mars_orbit = Orbit::new(
///     1.524 * AU,           // semi-major axis
///     0.0934,               // eccentricity
///     Degrees::new(1.85),   // inclination
///     Degrees::new(49.56),  // longitude of ascending node
///     Degrees::new(286.5),  // argument of perihelion
///     Degrees::new(19.41),  // mean anomaly at epoch
///     JulianDate::J2000,    // epoch
/// );
///
/// let mars_params = BodycentricParams::heliocentric(mars_orbit);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodycentricParams {
    /// The Keplerian orbital elements of the body.
    pub orbit: Orbit,
    /// Which standard center the orbit is defined relative to.
    pub orbit_center: OrbitReferenceCenter,
}

impl BodycentricParams {
    /// Creates parameters for a body-centered coordinate system.
    ///
    /// # Arguments
    ///
    /// - `orbit`: The Keplerian orbital elements of the body.
    /// - `orbit_center`: Which center the orbit is defined relative to.
    pub const fn new(orbit: Orbit, orbit_center: OrbitReferenceCenter) -> Self {
        Self {
            orbit,
            orbit_center,
        }
    }

    /// Creates parameters for a body orbiting the Sun (heliocentric orbit).
    ///
    /// This is the most common case for planets, asteroids, and comets.
    pub const fn heliocentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Heliocentric)
    }

    /// Creates parameters for a body orbiting Earth (geocentric orbit).
    ///
    /// Use this for artificial satellites, the Moon, etc.
    pub const fn geocentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Geocentric)
    }

    /// Creates parameters for a body orbiting the solar system barycenter.
    pub const fn barycentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Barycentric)
    }
}

impl Default for BodycentricParams {
    /// Returns default parameters with a circular 1 AU heliocentric orbit.
    ///
    /// Note: This default is primarily for internal use. In practice, you should
    /// always provide meaningful orbital elements for body-centric calculations.
    fn default() -> Self {
        use crate::astro::JulianDate;
        use qtty::AstronomicalUnits;

        Self {
            orbit: Orbit::new(
                AstronomicalUnits::new(1.0),
                0.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            ),
            orbit_center: OrbitReferenceCenter::Heliocentric,
        }
    }
}

/// Generic center for any orbiting celestial body.
///
/// This allows defining coordinate systems centered on satellites, planets,
/// moons, comets, asteroids, or any other body with known orbital elements.
/// The body's position at any time is computed via Keplerian propagation.
///
/// # Type Aliases
///
/// For convenience, you can use `Bodycentric` for any case, but conceptually:
/// - Satellites orbiting Earth: `Bodycentric` with `OrbitReferenceCenter::Geocentric`
/// - Planets/asteroids/comets: `Bodycentric` with `OrbitReferenceCenter::Heliocentric`
/// - Moons of other planets: Requires hierarchical orbit handling (future work)
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::algebra::centers::{Bodycentric, BodycentricParams, ReferenceCenter};
/// use siderust::coordinates::algebra::cartesian::Position;
/// use siderust::coordinates::algebra::frames;
/// use siderust::astro::orbit::Orbit;
/// use siderust::astro::JulianDate;
/// use qtty::*;
///
/// // Create orbital parameters for an Earth-orbiting satellite
/// let satellite_orbit = Orbit::new(
///     0.0000426 * AU,       // ~6378 km (low Earth orbit) in AU
///     0.001,                // nearly circular
///     Degrees::new(51.6),   // ISS-like inclination
///     Degrees::new(0.0),
///     Degrees::new(0.0),
///     Degrees::new(0.0),
///     JulianDate::J2000,
/// );
///
/// let sat_params = BodycentricParams::geocentric(satellite_orbit);
/// ```
#[derive(Debug, Copy, Clone)]
pub struct Bodycentric;

impl ReferenceCenter for Bodycentric {
    type Params = BodycentricParams;

    fn center_name() -> &'static str {
        "Bodycentric"
    }
}

// =============================================================================
// AffineCenter implementations for astronomical centers
// =============================================================================

// Implement AffineCenter for all actual coordinate centers.
// (AffineCenter trait itself comes from affn.)
impl AffineCenter for Barycentric {}
impl AffineCenter for Heliocentric {}
impl AffineCenter for Geocentric {}
impl AffineCenter for Topocentric {}
impl AffineCenter for Bodycentric {}


#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    #[test]
    fn center_names_are_correct() {
        assert_eq!(Barycentric::center_name(), "Barycentric");
        assert_eq!(Heliocentric::center_name(), "Heliocentric");
        assert_eq!(Topocentric::center_name(), "Topocentric");
        assert_eq!(Geocentric::center_name(), "Geocentric");
        assert_eq!(Bodycentric::center_name(), "Bodycentric");
        assert_eq!(<() as ReferenceCenter>::center_name(), "");
    }

    #[test]
    fn standard_centers_have_unit_params() {
        // Verify that standard centers use () as Params (zero-cost)
        let _: <Barycentric as ReferenceCenter>::Params = ();
        let _: <Heliocentric as ReferenceCenter>::Params = ();
        let _: <Geocentric as ReferenceCenter>::Params = ();
        let _: <() as ReferenceCenter>::Params = ();

        // Verify zero size
        assert_eq!(
            std::mem::size_of::<<Barycentric as ReferenceCenter>::Params>(),
            0
        );
        assert_eq!(
            std::mem::size_of::<<Heliocentric as ReferenceCenter>::Params>(),
            0
        );
        assert_eq!(
            std::mem::size_of::<<Geocentric as ReferenceCenter>::Params>(),
            0
        );
    }

    #[test]
    fn topocentric_has_observer_site_params() {
        // Verify Topocentric uses ObserverSite as Params
        let site = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let _: <Topocentric as ReferenceCenter>::Params = site;

        // Verify non-zero size (stores actual data)
        assert!(std::mem::size_of::<<Topocentric as ReferenceCenter>::Params>() > 0);
    }

    #[test]
    fn observer_site_default() {
        let site = ObserverSite::default();
        assert_eq!(site.lon.value(), 0.0);
        assert_eq!(site.lat.value(), 0.0);
        assert_eq!(site.height.value(), 0.0);
    }

    #[test]
    fn observer_site_equality() {
        let site1 = ObserverSite::new(10.0 * DEG, 20.0 * DEG, 100.0 * M);
        let site2 = ObserverSite::new(10.0 * DEG, 20.0 * DEG, 100.0 * M);
        let site3 = ObserverSite::new(10.0 * DEG, 20.0 * DEG, 200.0 * M);

        assert_eq!(site1, site2);
        assert_ne!(site1, site3);
    }

    #[test]
    fn bodycentric_has_params() {
        use crate::astro::JulianDate;

        // Create a simple orbit
        let orbit = Orbit::new(
            1.524 * AU,
            0.0934,
            Degrees::new(1.85),
            Degrees::new(49.56),
            Degrees::new(286.5),
            Degrees::new(19.41),
            JulianDate::J2000,
        );

        let params = BodycentricParams::heliocentric(orbit);
        let _: <Bodycentric as ReferenceCenter>::Params = params;

        // Verify non-zero size (stores actual data)
        assert!(std::mem::size_of::<<Bodycentric as ReferenceCenter>::Params>() > 0);
    }

    #[test]
    fn bodycentric_params_constructors() {
        use crate::astro::JulianDate;

        let orbit = Orbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let helio = BodycentricParams::heliocentric(orbit);
        assert_eq!(helio.orbit_center, OrbitReferenceCenter::Heliocentric);

        let geo = BodycentricParams::geocentric(orbit);
        assert_eq!(geo.orbit_center, OrbitReferenceCenter::Geocentric);

        let bary = BodycentricParams::barycentric(orbit);
        assert_eq!(bary.orbit_center, OrbitReferenceCenter::Barycentric);
    }

    #[test]
    fn bodycentric_params_default() {
        let params = BodycentricParams::default();
        assert_eq!(params.orbit_center, OrbitReferenceCenter::Heliocentric);
        assert_eq!(params.orbit.eccentricity, 0.0);
    }

    #[test]
    fn bodycentric_params_equality() {
        use crate::astro::JulianDate;

        let orbit1 = Orbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        let orbit2 = Orbit::new(
            2.0 * AU, // Different semi-major axis
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let params1 = BodycentricParams::heliocentric(orbit1);
        let params2 = BodycentricParams::heliocentric(orbit1);
        let params3 = BodycentricParams::heliocentric(orbit2);

        assert_eq!(params1, params2);
        assert_ne!(params1, params3);
    }
}
