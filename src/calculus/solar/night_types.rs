// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Night-related solar types (twilight thresholds, etc.).

use qtty::Degrees;

/// Common twilight types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Twilight {
    Civil,
    Nautical,
    Astronomical,
    Horizon,
    ApparentHorizon,
}

impl From<Twilight> for Degrees {
    fn from(t: Twilight) -> Degrees {
        match t {
            Twilight::Civil => Degrees::new(-6.0),
            Twilight::Nautical => Degrees::new(-12.0),
            Twilight::Astronomical => Degrees::new(-18.0),
            Twilight::Horizon => Degrees::new(0.0),
            Twilight::ApparentHorizon => Degrees::new(-0.833),
        }
    }
}

/// Standard twilight threshold definitions (Sun center altitude).
pub mod twilight {
    use qtty::Degrees;

    /// Civil twilight: Sun center 6° below horizon (-6°)
    pub const CIVIL: Degrees = Degrees::new(-6.0);

    /// Nautical twilight: Sun center 12° below horizon (-12°)
    pub const NAUTICAL: Degrees = Degrees::new(-12.0);

    /// Astronomical twilight: Sun center 18° below horizon (-18°)
    pub const ASTRONOMICAL: Degrees = Degrees::new(-18.0);

    /// Sunrise/sunset: Sun center at geometric horizon (0°)
    /// Note: For apparent sunrise/sunset, use -0.833° to account for refraction
    pub const HORIZON: Degrees = Degrees::new(0.0);

    /// Apparent sunrise/sunset accounting for atmospheric refraction (-0.833°)
    pub const APPARENT_HORIZON: Degrees = Degrees::new(-0.833);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_twilight_is_minus_6() {
        let d: Degrees = Twilight::Civil.into();
        assert!((d.value() - (-6.0)).abs() < 1e-12);
    }

    #[test]
    fn nautical_twilight_is_minus_12() {
        let d: Degrees = Twilight::Nautical.into();
        assert!((d.value() - (-12.0)).abs() < 1e-12);
    }

    #[test]
    fn astronomical_twilight_is_minus_18() {
        let d: Degrees = Twilight::Astronomical.into();
        assert!((d.value() - (-18.0)).abs() < 1e-12);
    }

    #[test]
    fn horizon_is_zero() {
        let d: Degrees = Twilight::Horizon.into();
        assert!((d.value() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn apparent_horizon_is_minus_0_833() {
        let d: Degrees = Twilight::ApparentHorizon.into();
        assert!((d.value() - (-0.833)).abs() < 1e-12);
    }

    #[test]
    fn twilight_clone_and_eq() {
        assert_eq!(Twilight::Civil, Twilight::Civil.clone());
        assert_ne!(Twilight::Civil, Twilight::Nautical);
    }

    #[test]
    fn twilight_debug() {
        let s = format!("{:?}", Twilight::Astronomical);
        assert!(s.contains("Astronomical"));
    }

    #[test]
    fn twilight_constants_match_variants() {
        assert!((twilight::CIVIL.value() - Degrees::from(Twilight::Civil).value()).abs() < 1e-12);
        assert!(
            (twilight::NAUTICAL.value() - Degrees::from(Twilight::Nautical).value()).abs() < 1e-12
        );
        assert!(
            (twilight::ASTRONOMICAL.value() - Degrees::from(Twilight::Astronomical).value()).abs()
                < 1e-12
        );
        assert!((twilight::HORIZON.value() - 0.0).abs() < 1e-12);
        assert!((twilight::APPARENT_HORIZON.value() - (-0.833)).abs() < 1e-12);
    }
}
