//! Frame-specific Direction constructors with astronomical conventions.

use super::direction_core::Direction;
use super::{clamp_polar, normalize_azimuth};
use crate::coordinates::frames;
use qtty::Degrees;

/// Macro to generate frame-specific Direction implementations with custom field names.
///
/// This macro reduces boilerplate by generating `new()`, and getter methods for different
/// astronomical reference frames, each with their own conventional field names.
///
/// # Syntax
/// ```ignore
/// impl_direction_frame! {
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
macro_rules! impl_direction_frame {
    ($frame:ty => ($az_name:ident, $polar_name:ident, $az_sym:expr, $polar_sym:expr, $desc:expr)) => {
        impl Direction<$frame> {
            #[doc = concat!("Creates a ", $desc, " direction.\n\n")]
            #[doc = "# Arguments\n"]
            #[doc = concat!("- `", stringify!($az_name), "`: ", $az_sym, ", normalized to [0°, 360°)\n")]
            #[doc = concat!("- `", stringify!($polar_name), "`: ", $polar_sym, ", clamped to [-90°, +90°]")]
            #[inline]
            pub fn new($az_name: Degrees, $polar_name: Degrees) -> Self {
                affn::spherical::Direction::new(
                    clamp_polar($polar_name),
                    normalize_azimuth($az_name),
                ).into()
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

impl_direction_frame!(
    frames::ICRS => (
        ra, dec,
        "Right Ascension (α), eastward from vernal equinox",
        "Declination (δ), from celestial equator",
        "ICRS"
    )
);

// =============================================================================
// Equatorial Frames
// =============================================================================

impl_direction_frame!(
    frames::EquatorialMeanJ2000 => (
        ra, dec,
        "Right Ascension (α), eastward from vernal equinox",
        "Declination (δ), from celestial equator",
        "mean-J2000 equatorial"
    )
);

impl_direction_frame!(
    frames::EquatorialMeanOfDate => (
        ra, dec,
        "Right Ascension (α), eastward from vernal equinox",
        "Declination (δ), from celestial equator",
        "mean-of-date equatorial"
    )
);

impl_direction_frame!(
    frames::EquatorialTrueOfDate => (
        ra, dec,
        "Right Ascension (α), eastward from vernal equinox",
        "Declination (δ), from celestial equator",
        "true-of-date equatorial"
    )
);

// =============================================================================
// Ecliptic
// =============================================================================

impl_direction_frame!(
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

// Custom implementation for Horizontal to follow IAU Alt-Az convention (altitude first)
impl Direction<frames::Horizontal> {
    /// Creates a Horizontal (Alt-Az) direction (IAU Alt-Az convention: altitude first).
    ///
    /// # Arguments
    /// - `alt`: Altitude (elevation) above the horizon, clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    #[inline]
    pub fn new(alt: Degrees, az: Degrees) -> Self {
        affn::spherical::Direction::new(clamp_polar(alt), normalize_azimuth(az)).into()
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

impl_direction_frame!(
    frames::ECEF => (
        lon, lat,
        "Geodetic longitude, positive eastward",
        "Geodetic latitude, positive northward",
        "Geographic (ECEF)"
    )
);
