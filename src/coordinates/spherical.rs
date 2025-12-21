//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! This module re-exports the algebraic types from `affn` and provides
//! frame-specific extension traits with convenient constructors and accessors
//! (e.g., `new_icrs(ra, dec)`, `ra()`, `dec()`).
//!
//! # Usage
//!
//! Import the extension traits to access frame-specific methods:
//!
//! ```rust
//! use siderust::coordinates::spherical::{direction, IcrsDirectionExt};
//! use qtty::*;
//!
//! // Extension trait method for ICRS
//! let dir = direction::ICRS::new_icrs(120.0 * DEG, 45.0 * DEG);
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! assert_eq!(dir.dec(), 45.0 * DEG);
//! ```

use crate::coordinates::{centers, frames};

// Re-export core types from affn
pub use affn::spherical::{Direction, Position};

use qtty::{Degrees, LengthUnit, Quantity};

// =============================================================================
// Direction type aliases (frame-only, no center)
// =============================================================================

pub mod direction {
    pub use super::Direction;
    use super::frames;

    /// **Ecliptic** direction (longitude *L*, latitude *B*).
    pub type Ecliptic = Direction<frames::Ecliptic>;
    /// **Equatorial** direction (right‑ascension *α*, declination *δ*).
    pub type Equatorial = Direction<frames::Equatorial>;
    /// **Horizontal** direction (azimuth *Az*, altitude *Alt*).
    pub type Horizontal = Direction<frames::Horizontal>;
    /// **ICRS** direction.
    pub type ICRS = Direction<frames::ICRS>;
    /// **Geographic** (ECEF) direction: latitude, longitude.
    pub type Geographic = Direction<frames::ECEF>;
}

// =============================================================================
// Position type aliases (center + frame + distance)
// =============================================================================

pub mod position {
    pub use super::Position;
    use super::{centers, frames};
    use qtty::Kilometer;

    /// **Heliocentric Ecliptic** coordinates *(L, B, R)*.
    ///
    /// * `L` – ecliptic longitude, degrees in `[0, 360)`
    /// * `B` – ecliptic latitude,  degrees in `[-90, 90]`
    /// * `R` – heliocentric distance in unit `U` (e.g. `AstronomicalUnit`)
    pub type Ecliptic<U> = Position<centers::Heliocentric, frames::Ecliptic, U>;

    /// **Geocentric Equatorial** coordinates *(δ, α, d)*.
    ///
    /// * `δ` – declination, degrees in `[-90, 90]`
    /// * `α` – right‑ascension, degrees in `[0, 360)`
    /// * `d` – geocentric distance in unit `U` (e.g. `Kilometer`)
    pub type Equatorial<U> = Position<centers::Geocentric, frames::Equatorial, U>;

    /// **Topocentric Horizontal** coordinates *(Alt, Az, d)*.
    ///
    /// * `Alt` – altitude above the horizon, degrees in `[-90, 90]`
    /// * `Az`  – azimuth from the north, degrees in `[0, 360)`
    /// * `d`   – straight‑line distance from the observer in unit `U`
    pub type Horizontal<U> = Position<centers::Topocentric, frames::Horizontal, U>;

    /// **Barycentric ICRS** coordinates.
    pub type ICRS<U> = Position<centers::Barycentric, frames::ICRS, U>;
    /// **Heliocentric ICRS** coordinates.
    pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
    /// **Geocentric ICRS** coordinates.
    pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;

    /// **Geographic (ECEF)** position: latitude, longitude, altitude (km).
    pub type Geographic = Position<centers::Geocentric, frames::ECEF, Kilometer>;
}

// =============================================================================
// Extension Traits
// =============================================================================

/// Macro to define extension traits for spherical coordinate systems.
///
/// This macro generates both `DirectionExt` and `PositionExt` traits for a given frame,
/// providing astronomy-specific constructors and accessors with appropriate naming.
///
/// # Parameters
///
/// - `$Frame`: The frame type (e.g., `ICRS`, `Equatorial`)
/// - `$prefix`: Trait name prefix (e.g., `Icrs` for `IcrsDirectionExt`)
/// - `$ctor`: Constructor method name (e.g., `new_icrs`)
/// - `$ctor_order`: Constructor argument order: either `azimuth_polar` or `polar_azimuth`
/// - `polar`: Mapping for the polar angle (name, doc string, field accessor)
/// - `azimuth`: Mapping for the azimuthal angle (name, doc string, field accessor)
macro_rules! define_spherical_extensions {
    (
        frame: $Frame:ty,
        prefix: $prefix:ident,
        constructor: $ctor:ident,
        ctor_order: azimuth_polar,
        polar: ($polar_name:ident, $polar_doc:expr, $polar_field:ident),
        azimuth: ($az_name:ident, $az_doc:expr, $az_field:ident)
        $(, special_position_impl: $special:tt)?
    ) => {
        paste::paste! {
            /// Extension trait for spherical directions.
            ///
            /// Provides astronomy-specific constructors and accessors.
            pub trait [<$prefix DirectionExt>] {
                #[doc = "Creates a new direction with " $az_doc " and " $polar_doc "."]
                fn $ctor($az_name: Degrees, $polar_name: Degrees) -> Self;
                
                #[doc = "Returns the " $az_doc " in degrees."]
                fn $az_name(&self) -> Degrees;
                
                #[doc = "Returns the " $polar_doc " in degrees."]
                fn $polar_name(&self) -> Degrees;
            }

            impl [<$prefix DirectionExt>] for Direction<$Frame> {
                #[inline]
                fn $ctor($az_name: Degrees, $polar_name: Degrees) -> Self {
                    Self::new($polar_name.wrap_quarter_fold(), $az_name.normalize())
                }
                
                #[inline]
                fn $az_name(&self) -> Degrees {
                    self.$az_field
                }
                
                #[inline]
                fn $polar_name(&self) -> Degrees {
                    self.$polar_field
                }
            }

            /// Extension trait for spherical positions.
            pub trait [<$prefix PositionExt>]<U: LengthUnit> {
                #[doc = "Creates a new position with " $az_doc ", " $polar_doc ", and distance."]
                fn $ctor<A: Into<Degrees>, T: Into<Quantity<U>>>($az_name: A, $polar_name: A, distance: T) -> Self;
                
                #[doc = "Returns the " $az_doc " in degrees."]
                fn $az_name(&self) -> Degrees;
                
                #[doc = "Returns the " $polar_doc " in degrees."]
                fn $polar_name(&self) -> Degrees;
            }

            impl<C: centers::ReferenceCenter<Params = ()>, U: LengthUnit> [<$prefix PositionExt>]<U>
                for Position<C, $Frame, U>
            {
                #[inline]
                fn $ctor<A: Into<Degrees>, T: Into<Quantity<U>>>($az_name: A, $polar_name: A, distance: T) -> Self {
                    Self::new_raw(
                        $polar_name.into().wrap_quarter_fold(),
                        $az_name.into().normalize(),
                        distance.into(),
                    )
                }
                
                #[inline]
                fn $az_name(&self) -> Degrees {
                    self.$az_field
                }
                
                #[inline]
                fn $polar_name(&self) -> Degrees {
                    self.$polar_field
                }
            }

            $(
                define_spherical_extensions!(@position_special $Frame, $prefix, $ctor, ($polar_name, $polar_doc, $polar_field), ($az_name, $az_doc, $az_field), $special);
            )?
        }
    };

    (
        frame: $Frame:ty,
        prefix: $prefix:ident,
        constructor: $ctor:ident,
        ctor_order: polar_azimuth,
        polar: ($polar_name:ident, $polar_doc:expr, $polar_field:ident),
        azimuth: ($az_name:ident, $az_doc:expr, $az_field:ident)
        $(, special_position_impl: $special:tt)?
    ) => {
        paste::paste! {
            /// Extension trait for spherical directions.
            ///
            /// Provides astronomy-specific constructors and accessors.
            pub trait [<$prefix DirectionExt>] {
                #[doc = "Creates a new direction with " $polar_doc " and " $az_doc "."]
                fn $ctor($polar_name: Degrees, $az_name: Degrees) -> Self;
                
                #[doc = "Returns the " $az_doc " in degrees."]
                fn $az_name(&self) -> Degrees;
                
                #[doc = "Returns the " $polar_doc " in degrees."]
                fn $polar_name(&self) -> Degrees;
            }

            impl [<$prefix DirectionExt>] for Direction<$Frame> {
                #[inline]
                fn $ctor($polar_name: Degrees, $az_name: Degrees) -> Self {
                    Self::new($polar_name.wrap_quarter_fold(), $az_name.normalize())
                }
                
                #[inline]
                fn $az_name(&self) -> Degrees {
                    self.$az_field
                }
                
                #[inline]
                fn $polar_name(&self) -> Degrees {
                    self.$polar_field
                }
            }

            /// Extension trait for spherical positions.
            pub trait [<$prefix PositionExt>]<U: LengthUnit> {
                #[doc = "Creates a new position with " $polar_doc ", " $az_doc ", and distance."]
                fn $ctor<A: Into<Degrees>, T: Into<Quantity<U>>>($polar_name: A, $az_name: A, distance: T) -> Self;
                
                #[doc = "Returns the " $az_doc " in degrees."]
                fn $az_name(&self) -> Degrees;
                
                #[doc = "Returns the " $polar_doc " in degrees."]
                fn $polar_name(&self) -> Degrees;
            }

            impl<C: centers::ReferenceCenter<Params = ()>, U: LengthUnit> [<$prefix PositionExt>]<U>
                for Position<C, $Frame, U>
            {
                #[inline]
                fn $ctor<A: Into<Degrees>, T: Into<Quantity<U>>>($polar_name: A, $az_name: A, distance: T) -> Self {
                    Self::new_raw(
                        $polar_name.into().wrap_quarter_fold(),
                        $az_name.into().normalize(),
                        distance.into(),
                    )
                }
                
                #[inline]
                fn $az_name(&self) -> Degrees {
                    self.$az_field
                }
                
                #[inline]
                fn $polar_name(&self) -> Degrees {
                    self.$polar_field
                }
            }

            $(
                define_spherical_extensions!(@position_special $Frame, $prefix, $ctor, ($polar_name, $polar_doc, $polar_field), ($az_name, $az_doc, $az_field), $special);
            )?
        }
    };

    // Special handling for position read-only traits (like HorizontalPositionReadExt)
    (@position_special $Frame:ty, $prefix:ident, $ctor:ident, 
     ($polar_name:ident, $polar_doc:expr, $polar_field:ident),
     ($az_name:ident, $az_doc:expr, $az_field:ident),
     read_topocentric) => {
        paste::paste! {
            /// Extension trait for reading position coordinates (works with Topocentric center).
            pub trait [<$prefix PositionReadExt>] {
                #[doc = "Returns the " $polar_doc " in degrees."]
                fn $polar_name(&self) -> Degrees;
                
                #[doc = "Returns the " $az_doc " in degrees."]
                fn $az_name(&self) -> Degrees;
            }

            impl<U: LengthUnit> [<$prefix PositionReadExt>]
                for Position<centers::Topocentric, $Frame, U>
            {
                #[inline]
                fn $polar_name(&self) -> Degrees {
                    self.$polar_field
                }
                
                #[inline]
                fn $az_name(&self) -> Degrees {
                    self.$az_field
                }
            }
        }
    };
}

// =============================================================================
// Extension Trait Definitions
// =============================================================================

define_spherical_extensions!(
    frame: frames::ICRS,
    prefix: Icrs,
    constructor: new_icrs,
    ctor_order: azimuth_polar,
    polar: (dec, "Declination (δ)", polar),
    azimuth: (ra, "Right Ascension (α)", azimuth)
);

define_spherical_extensions!(
    frame: frames::Equatorial,
    prefix: Equatorial,
    constructor: new_equatorial,
    ctor_order: azimuth_polar,
    polar: (dec, "Declination (δ)", polar),
    azimuth: (ra, "Right Ascension (α)", azimuth)
);

define_spherical_extensions!(
    frame: frames::Ecliptic,
    prefix: Ecliptic,
    constructor: new_ecliptic,
    ctor_order: azimuth_polar,
    polar: (lat, "ecliptic latitude (β)", polar),
    azimuth: (lon, "ecliptic longitude (λ)", azimuth)
);

define_spherical_extensions!(
    frame: frames::Horizontal,
    prefix: Horizontal,
    constructor: new_horizontal,
    ctor_order: polar_azimuth,
    polar: (alt, "altitude (elevation)", polar),
    azimuth: (az, "azimuth", azimuth),
    special_position_impl: read_topocentric
);

define_spherical_extensions!(
    frame: frames::ECEF,
    prefix: Ecef,
    constructor: new_geographic,
    ctor_order: polar_azimuth,
    polar: (lat, "geodetic latitude (φ)", polar),
    azimuth: (lon, "geodetic longitude (λ)", azimuth)
);
