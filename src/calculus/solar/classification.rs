// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Twilight Classification
//!
//! Classifies the current sky condition from the Sun's altitude angle, following
//! the IAU/USNO convention where each boundary belongs to the **higher** (brighter)
//! category — i.e. the upper bound is **inclusive** and the lower bound is
//! **exclusive** for every interval except the bottommost (`Dark`), which captures
//! all altitudes ≤ −18°.
//!
//! | Phase          | Condition                      |
//! |----------------|-------------------------------|
//! | `Day`          | sun_alt  >   0°               |
//! | `Civil`        | −6°  < sun_alt ≤   0°         |
//! | `Nautical`     | −12° < sun_alt ≤  −6°         |
//! | `Astronomical` | −18° < sun_alt ≤ −12°         |
//! | `Dark`         | sun_alt ≤ −18°                |
//!
//! ## Example
//!
//! ```rust
//! use siderust::qtty::Degrees;
//! use siderust::calculus::solar::classification::{TwilightPhase, twilight_classification};
//!
//! let phase = twilight_classification(Degrees::new(-7.5));
//! assert_eq!(phase, TwilightPhase::Nautical);
//! ```

use crate::qtty::{Angular, Deg, Quantity, Unit};

/// Sky condition derived from the Sun's altitude.
///
/// Boundaries are **inclusive on the upper (brighter) side**, following the
/// IAU/USNO convention: a Sun exactly at 0° is [`Day`](TwilightPhase::Day),
/// exactly at −6° is [`Civil`](TwilightPhase::Civil), etc.
///
/// See [`twilight_classification`] for the mapping function.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TwilightPhase {
    /// Sun above the horizon: altitude > 0°.
    Day,
    /// Civil twilight: −6° < altitude ≤ 0°.
    Civil,
    /// Nautical twilight: −12° < altitude ≤ −6°.
    Nautical,
    /// Astronomical twilight: −18° < altitude ≤ −12°.
    Astronomical,
    /// Full darkness: altitude ≤ −18°.
    Dark,
}

impl std::fmt::Display for TwilightPhase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Day => write!(f, "Day"),
            Self::Civil => write!(f, "Civil twilight"),
            Self::Nautical => write!(f, "Nautical twilight"),
            Self::Astronomical => write!(f, "Astronomical twilight"),
            Self::Dark => write!(f, "Dark"),
        }
    }
}

/// Classify the sky condition from the Sun's altitude.
///
/// The boundary convention follows IAU/USNO: the upper bound of each interval
/// is **inclusive** (`sun_alt == 0°` → [`Day`](TwilightPhase::Day);
/// `sun_alt == −6°` → [`Civil`](TwilightPhase::Civil); etc.).
///
/// ```text
/// sun_alt >   0°            → Day
/// −6° < sun_alt ≤   0°     → Civil
/// −12° < sun_alt ≤  −6°    → Nautical
/// −18° < sun_alt ≤ −12°    → Astronomical
/// sun_alt ≤ −18°            → Dark
/// ```
///
/// # Arguments
///
/// * `sun_altitude` — The Sun's apparent altitude as any typed angular
///   [`Quantity`] (e.g. [`Degrees`] or [`Radians`](crate::qtty::Radians)).
///   Negative values indicate the Sun is below the geometric horizon. The
///   value is converted to degrees internally, so callers may pass either
///   typed angle without manual conversion.
///
/// # Examples
///
/// ```rust
/// use siderust::qtty::{Degrees, Radians};
/// use siderust::calculus::solar::classification::{TwilightPhase, twilight_classification};
///
/// // Degrees:
/// assert_eq!(twilight_classification(Degrees::new(10.0)),  TwilightPhase::Day);
/// // 0° is the inclusive upper bound of Civil, so it classifies as Civil:
/// assert_eq!(twilight_classification(Degrees::new(0.0)),   TwilightPhase::Civil);
/// assert_eq!(twilight_classification(Degrees::new(-3.0)),  TwilightPhase::Civil);
/// // -6° is the inclusive upper bound of Nautical:
/// assert_eq!(twilight_classification(Degrees::new(-6.0)),  TwilightPhase::Nautical);
/// assert_eq!(twilight_classification(Degrees::new(-9.0)),  TwilightPhase::Nautical);
/// // -12° is the inclusive upper bound of Astronomical:
/// assert_eq!(twilight_classification(Degrees::new(-12.0)), TwilightPhase::Astronomical);
/// assert_eq!(twilight_classification(Degrees::new(-15.0)), TwilightPhase::Astronomical);
/// // -18° is the inclusive upper bound of Dark:
/// assert_eq!(twilight_classification(Degrees::new(-18.0)), TwilightPhase::Dark);
/// assert_eq!(twilight_classification(Degrees::new(-20.0)), TwilightPhase::Dark);
///
/// // Radians work the same — conversion is automatic:
/// let nautical_rad = Radians::new(-9.0_f64.to_radians());
/// assert_eq!(twilight_classification(nautical_rad), TwilightPhase::Nautical);
/// ```
pub fn twilight_classification<U>(sun_altitude: Quantity<U>) -> TwilightPhase
where
    U: Unit<Dim = Angular>,
{
    let alt = sun_altitude.to::<Deg>().value();
    if alt > 0.0 {
        TwilightPhase::Day
    } else if alt > -6.0 {
        TwilightPhase::Civil
    } else if alt > -12.0 {
        TwilightPhase::Nautical
    } else if alt > -18.0 {
        TwilightPhase::Astronomical
    } else {
        TwilightPhase::Dark
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Degrees;

    // ── Interior values ──────────────────────────────────────────────────────

    #[test]
    fn positive_altitude_is_day() {
        assert_eq!(
            twilight_classification(Degrees::new(0.001)),
            TwilightPhase::Day
        );
        assert_eq!(
            twilight_classification(Degrees::new(45.0)),
            TwilightPhase::Day
        );
        assert_eq!(
            twilight_classification(Degrees::new(90.0)),
            TwilightPhase::Day
        );
    }

    #[test]
    fn mid_civil_is_civil() {
        assert_eq!(
            twilight_classification(Degrees::new(-3.0)),
            TwilightPhase::Civil
        );
    }

    #[test]
    fn mid_nautical_is_nautical() {
        assert_eq!(
            twilight_classification(Degrees::new(-9.0)),
            TwilightPhase::Nautical
        );
    }

    #[test]
    fn mid_astronomical_is_astronomical() {
        assert_eq!(
            twilight_classification(Degrees::new(-15.0)),
            TwilightPhase::Astronomical
        );
    }

    #[test]
    fn deep_negative_is_dark() {
        assert_eq!(
            twilight_classification(Degrees::new(-30.0)),
            TwilightPhase::Dark
        );
        assert_eq!(
            twilight_classification(Degrees::new(-90.0)),
            TwilightPhase::Dark
        );
    }

    // ── Boundary values (upper-inclusive per IAU/USNO convention) ────────────

    #[test]
    fn boundary_zero_degrees_is_civil() {
        assert_eq!(
            twilight_classification(Degrees::new(0.0)),
            TwilightPhase::Civil,
            "0° is the inclusive upper bound of Civil (Day requires strictly > 0°)"
        );
    }

    #[test]
    fn just_above_zero_is_day() {
        assert_eq!(
            twilight_classification(Degrees::new(0.001)),
            TwilightPhase::Day
        );
    }

    #[test]
    fn just_below_zero_is_civil() {
        assert_eq!(
            twilight_classification(Degrees::new(-0.001)),
            TwilightPhase::Civil,
            "Just below 0° enters Civil"
        );
    }

    #[test]
    fn boundary_minus_6_is_nautical() {
        assert_eq!(
            twilight_classification(Degrees::new(-6.0)),
            TwilightPhase::Nautical,
            "-6° is the inclusive upper bound of Nautical"
        );
    }

    #[test]
    fn just_below_minus_6_is_nautical() {
        assert_eq!(
            twilight_classification(Degrees::new(-6.001)),
            TwilightPhase::Nautical,
            "Just below -6° enters Nautical"
        );
    }

    #[test]
    fn boundary_minus_12_is_astronomical() {
        assert_eq!(
            twilight_classification(Degrees::new(-12.0)),
            TwilightPhase::Astronomical,
            "-12° is the inclusive upper bound of Astronomical"
        );
    }

    #[test]
    fn just_below_minus_12_is_astronomical() {
        assert_eq!(
            twilight_classification(Degrees::new(-12.001)),
            TwilightPhase::Astronomical,
            "Just below -12° enters Astronomical"
        );
    }

    #[test]
    fn boundary_minus_18_is_dark() {
        assert_eq!(
            twilight_classification(Degrees::new(-18.0)),
            TwilightPhase::Dark,
            "-18° is the inclusive upper bound of Dark"
        );
    }

    #[test]
    fn just_below_minus_18_is_dark() {
        assert_eq!(
            twilight_classification(Degrees::new(-18.001)),
            TwilightPhase::Dark,
            "Just below -18° enters Dark"
        );
    }

    // ── Task-specified boundary cases ─────────────────────────────────────────

    #[test]
    fn task_boundary_zero() {
        // Task spec: Civil: -6° < sun_alt ≤ 0°  →  0° is Civil
        assert_eq!(
            twilight_classification(Degrees::new(0.0)),
            TwilightPhase::Civil
        );
    }

    #[test]
    fn task_boundary_minus_6() {
        // Task spec: Nautical: -12° < sun_alt ≤ -6°  →  -6° is Nautical
        assert_eq!(
            twilight_classification(Degrees::new(-6.0)),
            TwilightPhase::Nautical
        );
    }

    #[test]
    fn task_boundary_minus_12() {
        // Task spec: Astronomical: -18° < sun_alt ≤ -12°  →  -12° is Astronomical
        assert_eq!(
            twilight_classification(Degrees::new(-12.0)),
            TwilightPhase::Astronomical
        );
    }

    #[test]
    fn task_boundary_minus_18() {
        // Task spec: Dark: sun_alt ≤ -18°  →  -18° is Dark
        assert_eq!(
            twilight_classification(Degrees::new(-18.0)),
            TwilightPhase::Dark
        );
    }

    #[test]
    fn task_plus_0_001_is_day() {
        assert_eq!(
            twilight_classification(Degrees::new(0.001)),
            TwilightPhase::Day
        );
    }

    #[test]
    fn task_minus_18_001_is_dark() {
        assert_eq!(
            twilight_classification(Degrees::new(-18.001)),
            TwilightPhase::Dark
        );
    }

    #[test]
    fn radians_input_classifies_correctly() {
        use crate::qtty::Radians;
        // -9° → Nautical
        assert_eq!(
            twilight_classification(Radians::new(-9.0_f64.to_radians())),
            TwilightPhase::Nautical
        );
        // -18° boundary → Dark
        assert_eq!(
            twilight_classification(Radians::new(-18.0_f64.to_radians())),
            TwilightPhase::Dark
        );
        // +45° → Day
        assert_eq!(
            twilight_classification(Radians::new(45.0_f64.to_radians())),
            TwilightPhase::Day
        );
    }

    // ── Trait impls ───────────────────────────────────────────────────────────

    #[test]
    fn eq_and_clone() {
        assert_eq!(TwilightPhase::Civil, TwilightPhase::Civil.clone());
        assert_ne!(TwilightPhase::Civil, TwilightPhase::Nautical);
    }

    #[test]
    fn debug_contains_variant_name() {
        assert!(format!("{:?}", TwilightPhase::Dark).contains("Dark"));
        assert!(format!("{:?}", TwilightPhase::Astronomical).contains("Astronomical"));
    }

    #[test]
    fn display_is_human_readable() {
        assert_eq!(TwilightPhase::Day.to_string(), "Day");
        assert_eq!(TwilightPhase::Civil.to_string(), "Civil twilight");
        assert_eq!(TwilightPhase::Nautical.to_string(), "Nautical twilight");
        assert_eq!(
            TwilightPhase::Astronomical.to_string(),
            "Astronomical twilight"
        );
        assert_eq!(TwilightPhase::Dark.to_string(), "Dark");
    }

    #[test]
    fn hash_works_in_hashset() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(TwilightPhase::Day);
        set.insert(TwilightPhase::Civil);
        set.insert(TwilightPhase::Day); // duplicate
        assert_eq!(set.len(), 2);
    }
}
