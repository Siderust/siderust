// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed unit aliases for spacecraft dynamics quantities.
//!
//! ## Gravitational parameter
//!
//! The standard gravitational parameter `μ = G·M` is widely expressed in
//! **km³/s²** in astrodynamics (EGM2008 uses this convention).
//!
//! `GravitationalParameter` and the associated constants are now defined in
//! `qtty::dynamics` and re-exported here so that existing
//! `use crate::astro::dynamics::units::{GravitationalParameter, GM_EARTH, ...};`
//! paths continue to compile without change.

pub use crate::ext_qtty::dynamics::{
    GravitationalParameter, GravitationalParameterUnit, GM_EARTH, GM_MOON, GM_SUN,
};

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
