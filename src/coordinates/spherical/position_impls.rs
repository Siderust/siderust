//! Frame-specific Position constructors with astronomical conventions.

use super::{clamp_polar, normalize_azimuth};
use super::position_core::Position;
use crate::coordinates::{centers, frames};
use qtty::{Degrees, LengthUnit, Quantity};

/// Macro to generate frame-specific Position implementations with custom field names.
///
/// This macro reduces boilerplate by generating `new()` and getter methods for different
/// astronomical reference frames, each with their own conventional field names.
///
/// # Syntax
/// ```ignore
/// impl_position_frame! {
///     Frame => (azimuth_name, polar_name, azimuth_symbol, polar_symbol, description);
/// }
/// ```
///
/// # Arguments
/// - `Frame`: The reference frame type (e.g., `frames::ICRS`)
/// - `azimuth_name`: Name for the azimuth getter method (e.g., `ra`)
/// - `polar_name`: Name for the polar getter method (e.g., `dec`)
/// - `azimuth_symbol`: Symbol for documentation (e.g., `α`)
/// - `polar_symbol`: Symbol for documentation (e.g., `δ`)
/// - `description`: Description of the coordinate system
macro_rules! impl_position_frame {
    ($frame:ty => ($az_name:ident, $polar_name:ident, $az_sym:expr, $polar_sym:expr, $desc:expr)) => {
        impl<C, U> Position<C, $frame, U>
        where
            C: centers::ReferenceCenter<Params = ()>,
            U: LengthUnit,
        {
            #[doc = concat!("Creates a ", $desc, " position.\n\n")]
            #[doc = "# Arguments\n"]
            #[doc = concat!("- `", stringify!($az_name), "`: ", $az_sym, ", normalized to [0°, 360°)\n")]
            #[doc = concat!("- `", stringify!($polar_name), "`: ", $polar_sym, ", clamped to [-90°, +90°]\n")]
            #[doc = "- `distance`: Radial distance from the center"]
            #[inline]
            pub fn new<T: Into<Quantity<U>>>($az_name: Degrees, $polar_name: Degrees, distance: T) -> Self {
                Self::from_inner(affn::spherical::Position::new_raw(
                    clamp_polar($polar_name),
                    normalize_azimuth($az_name),
                    distance.into(),
                ))
            }

            #[doc = concat!("Returns the ", $az_sym, " in degrees.")]
            #[inline]
            pub fn $az_name(&self) -> Degrees {
                self.inner.azimuth
            }

            #[doc = concat!("Returns the ", $polar_sym, " in degrees.")]
            #[inline]
            pub fn $polar_name(&self) -> Degrees {
                self.inner.polar
            }
        }
    };
}

// =============================================================================
// ICRS
// =============================================================================

impl_position_frame!(
    frames::ICRS => (
        ra, dec,
        "Right Ascension (α)",
        "Declination (δ)",
        "ICRS"
    )
);

// =============================================================================
// Equatorial Frames
// =============================================================================

impl_position_frame!(
    frames::EquatorialMeanJ2000 => (
        ra, dec,
        "Right Ascension (α)",
        "Declination (δ)",
        "mean-J2000 equatorial"
    )
);

impl_position_frame!(
    frames::EquatorialMeanOfDate => (
        ra, dec,
        "Right Ascension (α)",
        "Declination (δ)",
        "mean-of-date equatorial"
    )
);

impl_position_frame!(
    frames::EquatorialTrueOfDate => (
        ra, dec,
        "Right Ascension (α)",
        "Declination (δ)",
        "true-of-date equatorial"
    )
);

// =============================================================================
// Ecliptic
// =============================================================================

impl_position_frame!(
    frames::Ecliptic => (
        lon, lat,
        "Ecliptic longitude (λ)",
        "Ecliptic latitude (β)",
        "Ecliptic"
    )
);

// =============================================================================
// Horizontal
// =============================================================================

// Note: Horizontal frame has reversed parameter order (alt, az) per IAU convention
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

impl_position_frame!(
    frames::ECEF => (
        lon, lat,
        "Geodetic longitude, positive eastward",
        "Geodetic latitude, positive northward",
        "Geographic (ECEF)"
    )
);
