//! # Reference Centers Module
//!
//! This module defines the concept of a *reference center* (origin) for astronomical and geodetic coordinate systems.
//! A reference center specifies the origin point from which positions are measured.
//!
//! ## Overview
//!
//! The [`ReferenceCenter`] trait provides a common interface for all reference center types. Each center is
//! represented as a zero-sized struct and implements the trait to provide its canonical name.
//!
//! The `new_center!` macro is used to conveniently declare new reference center types, ensuring consistency
//! and reducing boilerplate.
//!
//! ## Predefined Centers
//!
//! The following reference centers are provided out of the box:
//!
//! - `Barycentric`: Center of mass of the solar system.
//! - `Heliocentric`: Center of the Sun.
//! - `Geocentric`: Center of the Earth.
//! - `Topocentric`: Observer's location on the surface of the Earth (parameterized by [`ObserverSite`]).
//!
//! ## Parameterized Centers
//!
//! Some reference centers require runtime parameters. For example, [`Topocentric`] coordinates
//! need to know the observer's geographic location. This is achieved through the associated
//! `Params` type on [`ReferenceCenter`]:
//!
//! - For most centers (Barycentric, Heliocentric, Geocentric), `Params = ()` (zero-cost).
//! - For [`Topocentric`], `Params = ObserverSite` which stores longitude, latitude, and height.
//!
//! This allows coordinate values to carry their reference information inline without external context.
//!
//! ## Extending
//!
//! To define a new reference center, use the `new_center!` macro.
//!
//! This creates a new zero-sized type `Lunarcentric` that implements [`ReferenceCenter`].
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::centers::{ReferenceCenter, Geocentric};
//!
//! let name = Geocentric::center_name();
//! assert_eq!(name, "Geocentric");
//! ```

use qtty::{Degrees, Meter, Quantity};
use std::fmt::Debug;

/// A trait for defining a reference center (coordinate origin).
///
/// # Associated Types
///
/// - `Params`: Runtime parameters for this center. For most centers this is `()` (zero-cost).
///   For parameterized centers like [`Topocentric`], this carries observer location.
///
/// # Migration Note (v0.x → v0.y)
///
/// The `Params` associated type was added to support parameterized centers.
/// Existing code using `Barycentric`, `Heliocentric`, or `Geocentric` should continue
/// to work unchanged since their `Params = ()`.
pub trait ReferenceCenter {
    /// Runtime parameters for this center. Use `()` for centers that don't need parameters.
    type Params: Clone + Debug + Default + PartialEq;

    fn center_name() -> &'static str;
}

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
    /// use siderust::coordinates::centers::ObserverSite;
    /// use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
    ///
    /// let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    /// ```
    pub fn from_geographic(
        geo: &crate::coordinates::spherical::position::Geographic,
    ) -> Self {
        Self {
            lon: geo.azimuth,   // longitude is stored in the azimuth field
            lat: geo.polar,     // latitude is stored in the polar field
            height: geo.distance.to::<Meter>(),
        }
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
/// use siderust::coordinates::centers::{Topocentric, ObserverSite, ReferenceCenter};
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

impl ReferenceCenter for () {
    type Params = ();
    fn center_name() -> &'static str {
        ""
    }
}

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
}
