// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Horizontal Coordinate Convention Helpers
//!
//! This module provides explicit conversion utilities between common
//! horizontal (alt-az) coordinate conventions.
//!
//! ## Siderust's native convention
//!
//! Throughout siderust the **horizontal frame** uses:
//!
//! | Property        | Convention                              |
//! |-----------------|-----------------------------------------|
//! | **Azimuth origin** | North (0° = geographic north)         |
//! | **Azimuth sense**  | Clockwise (North → East → South → West) |
//! | **Altitude**       | Above the geometric horizon (positive upward) |
//!
//! This is the standard astronomical convention and matches IAU
//! recommendations.
//!
//! ## Other common conventions
//!
//! | Name                    | Origin | Sense            | Usage                                    |
//! |-------------------------|--------|------------------|------------------------------------------|
//! | **North-clockwise**     | North  | CW (N→E→S→W)    | Standard astronomy (siderust default)   |
//! | **South-clockwise**     | South  | CW (S→W→N→E)    | Some older European observatory logs     |
//! | **North-counter-clockwise** | North | CCW (N→W→S→E) | Mathematics / physics (ISO 31-11)        |
//! | **East-counter-clockwise**  | East  | CCW (E→N→W→S) | Meteorology, some engineering CAD tools  |
//!
//! ## Quick start
//!
//! ```rust
//! use siderust::coordinates::horizontal::{AzimuthOrigin, AzimuthSense, HorizontalConvention};
//! use siderust::coordinates::horizontal as hz;
//! use qtty::*;
//!
//! // Convert a foreign "south-clockwise" azimuth to siderust's native convention
//! let foreign_az = 45.0 * DEG; // 45° from South, clockwise
//! let native_az = hz::convert_azimuth(
//!     foreign_az,
//!     &HorizontalConvention::SOUTH_CLOCKWISE,
//!     &HorizontalConvention::NORTH_CLOCKWISE,
//! );
//! assert!((native_az.value() - 225.0).abs() < 1e-10);
//! ```
//!
//! The helpers work on raw [`Degrees`](qtty::Degrees) values as well as on
//! [`Direction<Horizontal>`] and [`Position<Topocentric, Horizontal, U>`]
//! coordinates.

use crate::coordinates::centers::Topocentric;
use crate::coordinates::frames::Horizontal;
use affn::spherical;
use qtty::{Degrees, LengthUnit, DEG};

// =============================================================================
// Convention descriptors
// =============================================================================

/// Cardinal point from which azimuth is counted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AzimuthOrigin {
    /// 0° = geographic north (siderust default).
    North,
    /// 0° = geographic south.
    South,
    /// 0° = geographic east.
    East,
    /// 0° = geographic west.
    West,
}

/// Direction in which azimuth increases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AzimuthSense {
    /// Azimuth grows clockwise when viewed from above (siderust default).
    Clockwise,
    /// Azimuth grows counter-clockwise when viewed from above.
    CounterClockwise,
}

/// A complete horizontal-frame azimuth convention.
///
/// Combines an [`AzimuthOrigin`] and an [`AzimuthSense`] so that a
/// conversion between any two conventions can be performed unambiguously.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HorizontalConvention {
    /// The cardinal point at which azimuth equals 0°.
    pub origin: AzimuthOrigin,
    /// The direction in which azimuth increases.
    pub sense: AzimuthSense,
}

impl HorizontalConvention {
    /// Create a new convention descriptor.
    pub const fn new(origin: AzimuthOrigin, sense: AzimuthSense) -> Self {
        Self { origin, sense }
    }

    /// Standard astronomical convention used by siderust.
    ///
    /// Azimuth measured from **North**, increasing **clockwise**
    /// (North → East → South → West).
    pub const NORTH_CLOCKWISE: Self = Self {
        origin: AzimuthOrigin::North,
        sense: AzimuthSense::Clockwise,
    };

    /// South-clockwise convention (some older European logs).
    ///
    /// Azimuth measured from **South**, increasing **clockwise**
    /// (South → West → North → East).
    pub const SOUTH_CLOCKWISE: Self = Self {
        origin: AzimuthOrigin::South,
        sense: AzimuthSense::Clockwise,
    };

    /// Mathematical / physics convention (ISO 31-11).
    ///
    /// Azimuth measured from **North**, increasing **counter-clockwise**
    /// (North → West → South → East).
    pub const NORTH_COUNTERCLOCKWISE: Self = Self {
        origin: AzimuthOrigin::North,
        sense: AzimuthSense::CounterClockwise,
    };

    /// Meteorological convention.
    ///
    /// Azimuth measured from **East**, increasing **counter-clockwise**
    /// (East → North → West → South).
    pub const EAST_COUNTERCLOCKWISE: Self = Self {
        origin: AzimuthOrigin::East,
        sense: AzimuthSense::CounterClockwise,
    };
}

impl std::fmt::Display for HorizontalConvention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let origin = match self.origin {
            AzimuthOrigin::North => "North",
            AzimuthOrigin::South => "South",
            AzimuthOrigin::East => "East",
            AzimuthOrigin::West => "West",
        };
        let sense = match self.sense {
            AzimuthSense::Clockwise => "clockwise",
            AzimuthSense::CounterClockwise => "counter-clockwise",
        };
        write!(f, "{origin}-{sense}")
    }
}

// =============================================================================
// Core azimuth conversion
// =============================================================================

/// Offset in degrees from North-clockwise for each origin, assuming CW sense.
///
/// This is the bearing of the origin cardinal point measured in the
/// standard North-CW system.
const fn origin_offset_cw(origin: AzimuthOrigin) -> f64 {
    match origin {
        AzimuthOrigin::North => 0.0,
        AzimuthOrigin::East => 90.0,
        AzimuthOrigin::South => 180.0,
        AzimuthOrigin::West => 270.0,
    }
}

/// Convert a raw azimuth value between two [`HorizontalConvention`]s.
///
/// The altitude component is unaffected by convention changes and is therefore
/// not part of this function.
///
/// # Algorithm
///
/// 1. Transform the input azimuth to an internal *North-clockwise*
///    representation by accounting for the source origin offset and sense.
/// 2. Transform from the internal representation to the target convention.
///
/// The result is always normalised to `[0°, 360°)`.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::horizontal::{HorizontalConvention, convert_azimuth};
/// use qtty::*;
///
/// // 90° in North-CW (= due East) → South-CW
/// let south_cw = convert_azimuth(
///     90.0 * DEG,
///     &HorizontalConvention::NORTH_CLOCKWISE,
///     &HorizontalConvention::SOUTH_CLOCKWISE,
/// );
/// assert!((south_cw.value() - 270.0).abs() < 1e-10);
/// ```
pub fn convert_azimuth(
    azimuth: Degrees,
    from: &HorizontalConvention,
    to: &HorizontalConvention,
) -> Degrees {
    if from == to {
        return azimuth.normalize();
    }

    // Step 1: source → North-CW
    let az_ncw = match from.sense {
        AzimuthSense::Clockwise => azimuth.value() + origin_offset_cw(from.origin),
        // CCW sense: negate first, then add origin offset
        AzimuthSense::CounterClockwise => -azimuth.value() + origin_offset_cw(from.origin),
    };

    // Step 2: North-CW → target
    let az_target = match to.sense {
        AzimuthSense::Clockwise => az_ncw - origin_offset_cw(to.origin),
        // CW → CCW: negate after subtracting target offset
        AzimuthSense::CounterClockwise => -(az_ncw - origin_offset_cw(to.origin)),
    };

    (az_target * DEG).normalize()
}

// =============================================================================
// Convenience: Direction<Horizontal> conversions
// =============================================================================

/// Convert a [`Direction<Horizontal>`] **from** a foreign convention
/// **to** siderust's native North-clockwise convention.
///
/// Altitude is preserved; only the azimuth is adjusted.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::horizontal::{HorizontalConvention, direction_to_native};
/// use siderust::coordinates::spherical::Direction;
/// use siderust::coordinates::frames::Horizontal;
/// use qtty::*;
///
/// let foreign = Direction::<Horizontal>::new(30.0 * DEG, 45.0 * DEG); // alt=30, az=45
/// let native = direction_to_native(&foreign, &HorizontalConvention::SOUTH_CLOCKWISE);
/// assert_eq!(native.alt(), 30.0 * DEG);
/// assert!((native.az().value() - 225.0).abs() < 1e-10);
/// ```
pub fn direction_to_native(
    dir: &spherical::Direction<Horizontal>,
    from: &HorizontalConvention,
) -> spherical::Direction<Horizontal> {
    let new_az = convert_azimuth(dir.az(), from, &HorizontalConvention::NORTH_CLOCKWISE);
    spherical::Direction::<Horizontal>::new(dir.alt(), new_az)
}

/// Convert a [`Direction<Horizontal>`] **from** siderust's native
/// North-clockwise convention **to** a foreign convention.
///
/// Altitude is preserved; only the azimuth is adjusted.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::horizontal::{HorizontalConvention, direction_from_native};
/// use siderust::coordinates::spherical::Direction;
/// use siderust::coordinates::frames::Horizontal;
/// use qtty::*;
///
/// let native = Direction::<Horizontal>::new(30.0 * DEG, 225.0 * DEG);
/// let foreign = direction_from_native(&native, &HorizontalConvention::SOUTH_CLOCKWISE);
/// assert_eq!(foreign.alt(), 30.0 * DEG);
/// assert!((foreign.az().value() - 45.0).abs() < 1e-10);
/// ```
pub fn direction_from_native(
    dir: &spherical::Direction<Horizontal>,
    to: &HorizontalConvention,
) -> spherical::Direction<Horizontal> {
    let new_az = convert_azimuth(dir.az(), &HorizontalConvention::NORTH_CLOCKWISE, to);
    spherical::Direction::<Horizontal>::new(dir.alt(), new_az)
}

/// Convert a [`Direction<Horizontal>`] between two arbitrary conventions.
///
/// Altitude is preserved; only the azimuth is adjusted.
pub fn convert_direction(
    dir: &spherical::Direction<Horizontal>,
    from: &HorizontalConvention,
    to: &HorizontalConvention,
) -> spherical::Direction<Horizontal> {
    let new_az = convert_azimuth(dir.az(), from, to);
    spherical::Direction::<Horizontal>::new(dir.alt(), new_az)
}

// =============================================================================
// Convenience: Position<Topocentric, Horizontal, U> conversions
// =============================================================================

/// Convert a [`Position<Topocentric, Horizontal, U>`] **from** a foreign
/// convention **to** siderust's native North-clockwise convention.
///
/// Altitude and distance are preserved; only the azimuth is adjusted.
/// The observer site (`center_params`) is also preserved.
pub fn position_to_native<U: LengthUnit>(
    pos: &spherical::Position<Topocentric, Horizontal, U>,
    from: &HorizontalConvention,
) -> spherical::Position<Topocentric, Horizontal, U> {
    let new_az = convert_azimuth(pos.az(), from, &HorizontalConvention::NORTH_CLOCKWISE);
    spherical::Position::<Topocentric, Horizontal, U>::new_raw_with_params(
        *pos.center_params(),
        pos.alt(),
        new_az,
        pos.distance,
    )
}

/// Convert a [`Position<Topocentric, Horizontal, U>`] **from**
/// siderust's native convention **to** a foreign convention.
///
/// Altitude and distance are preserved; only the azimuth is adjusted.
/// The observer site (`center_params`) is also preserved.
pub fn position_from_native<U: LengthUnit>(
    pos: &spherical::Position<Topocentric, Horizontal, U>,
    to: &HorizontalConvention,
) -> spherical::Position<Topocentric, Horizontal, U> {
    let new_az = convert_azimuth(pos.az(), &HorizontalConvention::NORTH_CLOCKWISE, to);
    spherical::Position::<Topocentric, Horizontal, U>::new_raw_with_params(
        *pos.center_params(),
        pos.alt(),
        new_az,
        pos.distance,
    )
}

/// Convert a [`Position<Topocentric, Horizontal, U>`] between two
/// arbitrary conventions.
///
/// Altitude and distance are preserved; only the azimuth is adjusted.
/// The observer site (`center_params`) is also preserved.
pub fn convert_position<U: LengthUnit>(
    pos: &spherical::Position<Topocentric, Horizontal, U>,
    from: &HorizontalConvention,
    to: &HorizontalConvention,
) -> spherical::Position<Topocentric, Horizontal, U> {
    let new_az = convert_azimuth(pos.az(), from, to);
    spherical::Position::<Topocentric, Horizontal, U>::new_raw_with_params(
        *pos.center_params(),
        pos.alt(),
        new_az,
        pos.distance,
    )
}

// =============================================================================
// Shorthand helpers for the most common "south ↔ north" swap
// =============================================================================

/// Flip an azimuth between North-origin and South-origin (both clockwise).
///
/// This is equivalent to adding/subtracting 180° and is the most common
/// convention mismatch encountered in practice.
///
/// The result is normalised to `[0°, 360°)`.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::horizontal::flip_north_south;
/// use qtty::*;
///
/// let north_az = 45.0 * DEG;   // NE in North-CW
/// let south_az = flip_north_south(north_az);
/// assert!((south_az.value() - 225.0).abs() < 1e-10);
///
/// // Round-trip
/// assert!((flip_north_south(south_az).value() - 45.0).abs() < 1e-10);
/// ```
pub fn flip_north_south(azimuth: Degrees) -> Degrees {
    (azimuth + 180.0 * DEG).normalize()
}

/// Reverse the sense of an azimuth (clockwise ↔ counter-clockwise)
/// while keeping the same origin.
///
/// This is equivalent to negating the angle (and re-normalising).
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::horizontal::flip_sense;
/// use qtty::*;
///
/// let cw_az = 90.0 * DEG;   // East in N-CW
/// let ccw_az = flip_sense(cw_az);
/// assert!((ccw_az.value() - 270.0).abs() < 1e-10);
/// ```
pub fn flip_sense(azimuth: Degrees) -> Degrees {
    (-(azimuth.value()) * DEG).normalize()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::DEG;

    const TOL: f64 = 1e-10;

    fn assert_az_eq(actual: Degrees, expected: f64) {
        let diff = (actual.value() - expected).abs();
        assert!(
            diff < TOL,
            "azimuth mismatch: got {}, expected {}",
            actual.value(),
            expected
        );
    }

    // -----------------------------------------------------------------
    // Identity
    // -----------------------------------------------------------------

    #[test]
    fn identity_north_cw() {
        let az = 123.456 * DEG;
        let result = convert_azimuth(
            az,
            &HorizontalConvention::NORTH_CLOCKWISE,
            &HorizontalConvention::NORTH_CLOCKWISE,
        );
        assert_az_eq(result, 123.456);
    }

    // -----------------------------------------------------------------
    // North-CW ↔ South-CW
    // -----------------------------------------------------------------

    #[test]
    fn north_cw_to_south_cw() {
        // 0° N-CW (North) → 180° S-CW
        assert_az_eq(
            convert_azimuth(
                0.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::SOUTH_CLOCKWISE,
            ),
            180.0,
        );
        // 90° N-CW (East) → 270° S-CW
        assert_az_eq(
            convert_azimuth(
                90.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::SOUTH_CLOCKWISE,
            ),
            270.0,
        );
        // 225° N-CW (SW) → 45° S-CW
        assert_az_eq(
            convert_azimuth(
                225.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::SOUTH_CLOCKWISE,
            ),
            45.0,
        );
    }

    #[test]
    fn south_cw_to_north_cw() {
        assert_az_eq(
            convert_azimuth(
                45.0 * DEG,
                &HorizontalConvention::SOUTH_CLOCKWISE,
                &HorizontalConvention::NORTH_CLOCKWISE,
            ),
            225.0,
        );
    }

    // -----------------------------------------------------------------
    // North-CW ↔ North-CCW
    // -----------------------------------------------------------------

    #[test]
    fn north_cw_to_north_ccw() {
        // 90° N-CW (East) → 270° N-CCW
        assert_az_eq(
            convert_azimuth(
                90.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::NORTH_COUNTERCLOCKWISE,
            ),
            270.0,
        );
        // 0° N-CW → 0° N-CCW (North stays north)
        assert_az_eq(
            convert_azimuth(
                0.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::NORTH_COUNTERCLOCKWISE,
            ),
            0.0,
        );
    }

    #[test]
    fn north_ccw_to_north_cw() {
        // 90° N-CCW (West) → 270° N-CW
        assert_az_eq(
            convert_azimuth(
                90.0 * DEG,
                &HorizontalConvention::NORTH_COUNTERCLOCKWISE,
                &HorizontalConvention::NORTH_CLOCKWISE,
            ),
            270.0,
        );
    }

    // -----------------------------------------------------------------
    // North-CW ↔ East-CCW (meteorological)
    // -----------------------------------------------------------------

    #[test]
    fn north_cw_to_east_ccw() {
        // 0° N-CW (North) → 90° E-CCW
        assert_az_eq(
            convert_azimuth(
                0.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::EAST_COUNTERCLOCKWISE,
            ),
            90.0,
        );
        // 90° N-CW (East) → 0° E-CCW
        assert_az_eq(
            convert_azimuth(
                90.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::EAST_COUNTERCLOCKWISE,
            ),
            0.0,
        );
        // 180° N-CW (South) → 270° E-CCW
        assert_az_eq(
            convert_azimuth(
                180.0 * DEG,
                &HorizontalConvention::NORTH_CLOCKWISE,
                &HorizontalConvention::EAST_COUNTERCLOCKWISE,
            ),
            270.0,
        );
    }

    #[test]
    fn east_ccw_to_north_cw() {
        // 0° E-CCW (East) → 90° N-CW
        assert_az_eq(
            convert_azimuth(
                0.0 * DEG,
                &HorizontalConvention::EAST_COUNTERCLOCKWISE,
                &HorizontalConvention::NORTH_CLOCKWISE,
            ),
            90.0,
        );
    }

    // -----------------------------------------------------------------
    // Round-trip invariant
    // -----------------------------------------------------------------

    #[test]
    fn roundtrip_all_conventions() {
        let conventions = [
            HorizontalConvention::NORTH_CLOCKWISE,
            HorizontalConvention::SOUTH_CLOCKWISE,
            HorizontalConvention::NORTH_COUNTERCLOCKWISE,
            HorizontalConvention::EAST_COUNTERCLOCKWISE,
        ];
        let test_azimuths = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 359.9];

        for &from in &conventions {
            for &to in &conventions {
                for &az_val in &test_azimuths {
                    let az = az_val * DEG;
                    let converted = convert_azimuth(az, &from, &to);
                    let back = convert_azimuth(converted, &to, &from);
                    assert!(
                        (back.value() - az.normalize().value()).abs() < TOL,
                        "Round-trip failed: {az_val}° from {from} → {to} → {from}: got {}",
                        back.value()
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // Shorthand helpers
    // -----------------------------------------------------------------

    #[test]
    fn flip_north_south_basic() {
        assert_az_eq(flip_north_south(0.0 * DEG), 180.0);
        assert_az_eq(flip_north_south(45.0 * DEG), 225.0);
        assert_az_eq(flip_north_south(180.0 * DEG), 0.0);
        assert_az_eq(flip_north_south(270.0 * DEG), 90.0);
    }

    #[test]
    fn flip_sense_basic() {
        assert_az_eq(flip_sense(90.0 * DEG), 270.0);
        assert_az_eq(flip_sense(0.0 * DEG), 0.0);
        assert_az_eq(flip_sense(45.0 * DEG), 315.0);
    }

    // -----------------------------------------------------------------
    // Direction helpers
    // -----------------------------------------------------------------

    #[test]
    fn direction_to_native_south_cw() {
        // A direction with alt=30°, az=45° in South-CW
        let foreign = spherical::Direction::<Horizontal>::new(30.0 * DEG, 45.0 * DEG);
        let native = direction_to_native(&foreign, &HorizontalConvention::SOUTH_CLOCKWISE);
        assert_eq!(native.alt(), 30.0 * DEG);
        assert_az_eq(native.az(), 225.0);
    }

    #[test]
    fn direction_from_native_south_cw() {
        let native = spherical::Direction::<Horizontal>::new(30.0 * DEG, 225.0 * DEG);
        let foreign = direction_from_native(&native, &HorizontalConvention::SOUTH_CLOCKWISE);
        assert_eq!(foreign.alt(), 30.0 * DEG);
        assert_az_eq(foreign.az(), 45.0);
    }

    #[test]
    fn direction_roundtrip() {
        let original = spherical::Direction::<Horizontal>::new(60.0 * DEG, 123.0 * DEG);
        let converted =
            convert_direction(&original, &HorizontalConvention::NORTH_CLOCKWISE, &HorizontalConvention::EAST_COUNTERCLOCKWISE);
        let back = convert_direction(
            &converted,
            &HorizontalConvention::EAST_COUNTERCLOCKWISE,
            &HorizontalConvention::NORTH_CLOCKWISE,
        );
        assert_eq!(back.alt(), original.alt());
        assert_az_eq(back.az(), original.az().value());
    }

    // -----------------------------------------------------------------
    // Display
    // -----------------------------------------------------------------

    #[test]
    fn convention_display() {
        assert_eq!(
            HorizontalConvention::NORTH_CLOCKWISE.to_string(),
            "North-clockwise"
        );
        assert_eq!(
            HorizontalConvention::EAST_COUNTERCLOCKWISE.to_string(),
            "East-counter-clockwise"
        );
    }
}
