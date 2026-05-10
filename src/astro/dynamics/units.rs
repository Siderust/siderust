// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed unit aliases for spacecraft dynamics quantities.
//!
//! ## Gravitational parameter
//!
//! The standard gravitational parameter `μ = G·M` is widely expressed in
//! **km³/s²** in astrodynamics (EGM2008 uses this convention). `qtty` does not
//! yet carry a named unit for this dimension; the alias below is defined
//! locally as `Per<CubicKilometer, Prod<Second, Second>>`.
//!
//! When `qtty` gains a `GravitationalParameter` unit this alias should be
//! replaced by the upstream type and this file may become empty.

use crate::qtty::unit::Second;
use crate::qtty::{CubicKilometer, Per, Prod, Quantity};

/// Typed gravitational parameter: km³ / s².
///
/// Used for `GM` of Earth, Sun, Moon, and other bodies throughout the
/// dynamics module. The numerical value matches the EGM2008 / WGS-84
/// convention (e.g. Earth GM ≈ 398 600.4418 km³/s²).
///
/// # Notes
///
/// `CubicKilometer` = `Prod<Prod<Kilometer, Kilometer>, Kilometer>` and
/// `Prod<Second, Second>` is the squared-second unit, so this is exactly
/// the `L³ T⁻²` dimension expected for a gravitational parameter.
pub type GravitationalParameter = Quantity<Per<CubicKilometer, Prod<Second, Second>>>;

/// Earth standard gravitational parameter (EGM2008 / WGS-84).
pub const GM_EARTH: GravitationalParameter = GravitationalParameter::new(398_600.441_8);

/// Sun standard gravitational parameter (IAU 2012 system).
pub const GM_SUN: GravitationalParameter = GravitationalParameter::new(1.327_124_400_18e11);

/// Moon standard gravitational parameter (DE430 value).
pub const GM_MOON: GravitationalParameter = GravitationalParameter::new(4.902_800_066e3);


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gm_earth_value() {
        assert!((GM_EARTH.value() - 398_600.441_8).abs() < 1e-3);
    }

    #[test]
    fn gm_sun_order_of_magnitude() {
        assert!(GM_SUN.value() > 1e11 && GM_SUN.value() < 1.4e11);
    }

    #[test]
    fn gm_moon_order_of_magnitude() {
        assert!(GM_MOON.value() > 4_000.0 && GM_MOON.value() < 6_000.0);
    }
}
