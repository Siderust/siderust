// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ICAO International Standard Atmosphere (ISA)
//!
//! ## Scientific scope
//!
//! The ICAO Standard Atmosphere (ISA) is the reference pressure-altitude
//! model defined in ICAO Doc 7488, 3rd edition (1993). It divides the
//! atmosphere into layers with constant temperature lapse rates; within each
//! layer pressure and density follow hydrostatic equilibrium.
//!
//! This module implements:
//! - Pressure from geometric altitude, covering the troposphere (0–11 km)
//!   and lower stratosphere (11–20 km).
//! - The inverse: pressure altitude from static pressure.
//! - Geopotential ↔ geometric altitude conversion using the standard WGS-84
//!   Earth radius of 6 356.766 km.
//!
//! Validity range: geometric altitude 0–20 km (−2 000 m to 20 000 m is
//! supported; values outside produce extrapolated results).
//!
//! ## Technical scope
//!
//! All public functions accept and return typed [`crate::qtty`] values:
//!
//! | Function | Input | Output |
//! |----------|-------|--------|
//! | [`pressure_at_altitude`] | geometric altitude ([`Meters`]) | static pressure ([`Pascals`]) |
//! | [`pressure_altitude_m`] | static pressure ([`Pascals`]) | pressure altitude ([`Meters`]) |
//! | [`temperature_at_altitude`] | geometric altitude ([`Meters`]) | temperature ([`Kelvins`]) |
//! | [`geopotential_altitude_m`] | geometric altitude ([`Meters`]) | geopotential altitude ([`Meters`]) |
//! | [`geometric_altitude_m`] | geopotential altitude ([`Meters`]) | geometric altitude ([`Meters`]) |
//!
//! ## References
//!
//! - ICAO (1993). *Manual of the ICAO Standard Atmosphere*, 3rd ed. Doc
//!   7488. International Civil Aviation Organization.
//! - ICAO (2010). *Annex 2 — Rules of the Air*, 10th edition.

use crate::qtty::{Kelvins, Meters, Pascals};

// =============================================================================
// ISA constants
// =============================================================================

/// ISA sea-level pressure (Pa).
const P0: f64 = 101_325.0;

/// ISA sea-level temperature (K).
const T0: f64 = 288.15;

/// ISA tropopause temperature (K); constant from 11 km to 20 km.
const T_TROPO: f64 = 216.65;

/// Tropospheric lapse rate (K/m).
const LAPSE_RATE: f64 = -6.5e-3;

/// Tropopause geometric altitude (m).
const H_TROPO: f64 = 11_000.0;

/// ISA molar mass of dry air (kg/mol).
const M: f64 = 0.028_964_4;

/// Universal gas constant (J/(mol·K)).
const R: f64 = 8.314_462_618;

/// Standard acceleration of gravity (m/s²).
const G0: f64 = 9.806_65;

/// WGS-84 mean Earth radius used for geopotential altitude (km → m: 6 356 766 m).
const R_EARTH: f64 = 6_356_766.0;

// Derived constant for the barometric formula: g0·M / (R·|L|).
const EXP_TROP: f64 = G0 * M / (R * (-LAPSE_RATE));

/// ISA pressure at the tropopause base (Pa).
///
/// Computed from the troposphere formula at H_TROPO to ensure consistency
/// with `pressure_at_altitude` and its inverse.
#[inline]
fn p_tropo() -> f64 {
    let t_ratio = 1.0 + LAPSE_RATE * H_TROPO / T0;
    P0 * t_ratio.powf(EXP_TROP)
}

// =============================================================================
// Public API
// =============================================================================

/// Static pressure at a given geometric altitude.
///
/// Troposphere (0–11 km): hydrostatic gradient model.
/// Lower stratosphere (11–20 km): isothermal at 216.65 K.
///
/// # Arguments
///
/// - `altitude` — geometric altitude above MSL.
///
/// # Returns
///
/// Static pressure in [`Pascals`].
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::aircraft::isa::pressure_at_altitude;
/// use siderust::qtty::Meters;
///
/// // Sea level: 101 325 Pa
/// let p_sl = pressure_at_altitude(Meters::new(0.0));
/// assert!((p_sl.value() - 101_325.0).abs() < 1.0);
///
/// // FL350 ≈ 10 668 m: ~23 842 Pa
/// let p_fl350 = pressure_at_altitude(Meters::new(10_668.0));
/// assert!(p_fl350.value() > 23_000.0 && p_fl350.value() < 24_500.0);
/// ```
pub fn pressure_at_altitude(altitude: Meters) -> Pascals {
    let h = altitude.value();
    let p = if h <= H_TROPO {
        // Troposphere: p = p0 * (T / T0) ^ (g0·M / (R·|L|))
        let t_ratio = 1.0 + LAPSE_RATE * h / T0;
        P0 * t_ratio.powf(EXP_TROP)
    } else {
        // Isothermal stratosphere: p = p_tropo * exp(−g0·M·Δh / (R·T_tropo))
        let delta_h = h - H_TROPO;
        p_tropo() * (-G0 * M * delta_h / (R * T_TROPO)).exp()
    };
    Pascals::new(p)
}

/// Pressure altitude (ISA barometric altitude) from static pressure.
///
/// Inverts [`pressure_at_altitude`].  Returns the geometric altitude at which
/// the ISA model produces the given pressure.
///
/// # Arguments
///
/// - `pressure` — static pressure in [`Pascals`].
///
/// # Returns
///
/// Pressure altitude in [`Meters`].
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::aircraft::isa::{pressure_at_altitude, pressure_altitude_m};
/// use siderust::qtty::Meters;
///
/// let h_in = Meters::new(8_500.0);
/// let p = pressure_at_altitude(h_in);
/// let h_out = pressure_altitude_m(p);
/// assert!((h_out.value() - h_in.value()).abs() < 0.01);
/// ```
pub fn pressure_altitude_m(pressure: Pascals) -> Meters {
    let p = pressure.value();
    let pt = p_tropo();
    let h = if p >= pt {
        // Troposphere inverse:
        //   h = (T0 / L) * ((p / p0)^(1 / EXP_TROP) − 1)
        //   L is negative so T0/L < 0, and (p/P0)^(1/EXP_TROP) − 1 < 0, giving h > 0.
        (T0 / LAPSE_RATE) * ((p / P0).powf(1.0 / EXP_TROP) - 1.0)
    } else {
        // Stratosphere inverse:
        //   h = H_TROPO − (R·T_tropo / (g0·M)) * ln(p / p_tropo)
        H_TROPO - (R * T_TROPO / (G0 * M)) * (p / pt).ln()
    };
    Meters::new(h)
}

/// ISA air temperature at a given geometric altitude.
///
/// # Arguments
///
/// - `altitude` — geometric altitude above MSL.
///
/// # Returns
///
/// Temperature in [`Kelvins`].
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::aircraft::isa::temperature_at_altitude;
/// use siderust::qtty::Meters;
///
/// // Sea level: 288.15 K
/// let t = temperature_at_altitude(Meters::new(0.0));
/// assert!((t.value() - 288.15).abs() < 1e-6);
///
/// // Above tropopause: isothermal at 216.65 K
/// let t_strat = temperature_at_altitude(Meters::new(15_000.0));
/// assert!((t_strat.value() - 216.65).abs() < 1e-6);
/// ```
pub fn temperature_at_altitude(altitude: Meters) -> Kelvins {
    let h = altitude.value();
    let t = if h <= H_TROPO {
        T0 + LAPSE_RATE * h
    } else {
        T_TROPO
    };
    Kelvins::new(t)
}

/// Geopotential altitude from geometric altitude.
///
/// Uses the standard WGS-84 mean Earth radius (6 356 766 m).
///
/// # Arguments
///
/// - `geometric` — geometric altitude above MSL.
///
/// # Returns
///
/// Geopotential altitude in [`Meters`].
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::aircraft::isa::geopotential_altitude_m;
/// use siderust::qtty::Meters;
///
/// // At sea level, geopotential ≈ geometric
/// let h = geopotential_altitude_m(Meters::new(0.0));
/// assert!((h.value() - 0.0).abs() < 1e-6);
///
/// // At 10 km the difference is ~15 m
/// let h = geopotential_altitude_m(Meters::new(10_000.0));
/// assert!((h.value() - 10_000.0).abs() < 20.0);
/// ```
pub fn geopotential_altitude_m(geometric: Meters) -> Meters {
    let z = geometric.value();
    Meters::new(R_EARTH * z / (R_EARTH + z))
}

/// Geometric altitude from geopotential altitude.
///
/// Inverse of [`geopotential_altitude_m`].
///
/// # Arguments
///
/// - `geopotential` — geopotential altitude above MSL.
///
/// # Returns
///
/// Geometric altitude in [`Meters`].
///
/// # Examples
///
/// ```rust
/// use siderust::bodies::aircraft::isa::{geopotential_altitude_m, geometric_altitude_m};
/// use siderust::qtty::Meters;
///
/// let z_in = Meters::new(10_000.0);
/// let hp = geopotential_altitude_m(z_in);
/// let z_out = geometric_altitude_m(hp);
/// assert!((z_out.value() - z_in.value()).abs() < 1e-6);
/// ```
pub fn geometric_altitude_m(geopotential: Meters) -> Meters {
    let hp = geopotential.value();
    Meters::new(R_EARTH * hp / (R_EARTH - hp))
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_pressure() {
        let p = pressure_at_altitude(Meters::new(0.0));
        assert!((p.value() - 101_325.0).abs() < 1.0, "p = {}", p.value());
    }

    #[test]
    fn tropopause_pressure() {
        let p = pressure_at_altitude(Meters::new(H_TROPO));
        let expected = p_tropo();
        assert!((p.value() - expected).abs() < 0.1, "p = {}", p.value());
    }

    #[test]
    fn stratosphere_pressure_20km() {
        let p = pressure_at_altitude(Meters::new(20_000.0));
        // ICAO Doc 7488: 5 474.88 Pa at 20 km
        assert!((p.value() - 5_474.88).abs() < 5.0, "p = {}", p.value());
    }

    #[test]
    fn pressure_altitude_roundtrip() {
        for h in [0.0_f64, 5_000.0, 11_000.0, 15_000.0, 20_000.0] {
            let p = pressure_at_altitude(Meters::new(h));
            let h_out = pressure_altitude_m(p);
            assert!(
                (h_out.value() - h).abs() < 0.01,
                "h_in={h}, h_out={}",
                h_out.value()
            );
        }
    }

    #[test]
    fn sea_level_temperature() {
        let t = temperature_at_altitude(Meters::new(0.0));
        assert!((t.value() - 288.15).abs() < 1e-6);
    }

    #[test]
    fn stratosphere_temperature_isothermal() {
        let t = temperature_at_altitude(Meters::new(15_000.0));
        assert!((t.value() - 216.65).abs() < 1e-6);
    }

    #[test]
    fn geopotential_geometric_roundtrip() {
        for h in [0.0_f64, 1_000.0, 10_000.0, 20_000.0] {
            let hp = geopotential_altitude_m(Meters::new(h));
            let z = geometric_altitude_m(hp);
            assert!((z.value() - h).abs() < 1e-6, "h={h}, z={}", z.value());
        }
    }
}
