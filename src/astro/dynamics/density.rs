// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Atmospheric density models for satellite drag.
//!
//! ## Scientific scope
//!
//! Satellite drag acceleration requires the atmospheric mass density `ρ`
//! (kg/m³) at the spacecraft's geodetic altitude.  This module provides a
//! trait-based abstraction over density sources so that force-model
//! implementations remain independent of the atmosphere backend.
//!
//! Two built-in models are provided:
//!
//! | Type | Accuracy | When to use |
//! |------|----------|-------------|
//! | [`ExponentialAtmosphere`] | ±factor-2 near reference altitude | Unit tests, quick closed-form propagation |
//! | [`Nrlmsise00LiteApprox`] | ±50 % at mean F10.7 = 140 | Teaching / regression baselines |
//! | [`ConstantDensity`] | Synthetic only | Controlled unit tests |
//!
//! Neither built-in model is suitable for operational precision orbit determination.
//! Plug in a full NRLMSISE-00 / JB2008 / DTM2000 backend by implementing
//! [`DensityProvider`] in a downstream crate.
//!
//! ## Geodetic altitude
//!
//! [`geodetic_altitude`] converts an ECEF (or GCRS — see caveat below)
//! Cartesian position to geodetic altitude above the WGS-84 ellipsoid using
//! the iterative Bowring (1985) method.  The current implementation uses a
//! **placeholder identity GCRS↔ITRF rotation** for drag altitude.  The
//! resulting altitude error is at most ~20 km at the poles, which is within
//! the uncertainty of every built-in density model.  The
//! [`Geopotential`][crate::astro::dynamics::forces::Geopotential] force model
//! applies the Earth Rotation Angle for non-zonal terms, but drag does not
//! yet carry an epoch-dependent rotation.
//!
//! ## Unit convention
//!
//! All density inputs/outputs use SI units (kg/m³).  The drag force model
//! works internally in km/s and converts:
//! `a [km/s²] = −½·Cd·(A/m)·ρ·|v_rel_SI|·v_rel`
//! where `ρ` in kg/m³, `A/m` in m²/kg, `|v_rel_SI|` in m/s, `v_rel` in km/s.
//!
//! ## References
//!
//! * Picone, J.M. et al. (2002), "NRLMSISE-00 empirical model of the atmosphere",
//!   *J. Geophys. Res.*, 107(A12), 1468.
//! * Vallado, D.A. (2013), *Fundamentals of Astrodynamics and Applications*, §8.6, §3.4.
//! * Montenbruck & Gill (2000), *Satellite Orbits*, §3.5.
//! * Bowring, B.R. (1985), "The geodetic line and the geodetic coordinates",
//!   *Survey Review*, 28, 276–281.

use affn::cartesian::Position;
use siderust_archive::atmosphere::tables::NRLMSISE_TABLE;

use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::unit::Kilometer;
use crate::qtty::{KilogramsPerCubicMeter, Kilometers, Ratios};

use super::errors::DynamicsError;

// =============================================================================
// DensityProvider trait
// =============================================================================

/// Provider returning atmospheric mass density `ρ` (kg/m³) at a given geodetic altitude.
///
/// Implement this trait to plug in any density model into the drag force model.
/// All implementing types must be [`Send`] + [`Sync`] so they can be shared
/// across threads via [`std::sync::Arc`].
///
/// # Error handling
///
/// Returns [`DynamicsError::AltitudeBelowSurface`] when `altitude < 0 km`.
/// Other implementation-specific failures should use
/// [`DynamicsError::AtmosphereProviderError`].
pub trait DensityProvider: Send + Sync {
    /// Human-readable name of this density model.
    fn name(&self) -> &'static str;

    /// Mass density `ρ` (kg/m³) at the given geodetic `altitude` (km).
    ///
    /// # Errors
    ///
    /// Returns [`DynamicsError::AltitudeBelowSurface`] if `altitude < 0 km`.
    fn density(&self, altitude: Kilometers) -> Result<KilogramsPerCubicMeter, DynamicsError>;
}

/// Convenience alias for dynamic [`DensityProvider`] trait objects.
///
/// Use `Arc<AtmosphereProvider>` to hold a heap-allocated, thread-safe
/// density provider in a [`DynamicsContext`].
///
/// [`DynamicsContext`]: super::context::DynamicsContext
pub type AtmosphereProvider = dyn DensityProvider + Send + Sync;

// =============================================================================
// ExponentialAtmosphere
// =============================================================================

/// Single-layer exponential atmosphere.
///
/// ```text
/// ρ(h) = ρ₀ · exp(−(h − h₀) / H)
/// ```
///
/// where:
/// - `ρ₀` is the reference density at reference altitude `h₀`,
/// - `H` is the scale height.
///
/// ## Accuracy
///
/// This profile is representative only within ~±100 km of the reference
/// altitude `h₀`; outside that band errors can exceed one order of magnitude.
/// Below 100 km the exponential extrapolation is physically meaningless.
///
/// ## References
///
/// * Vallado §8.6, Table 8-4.
#[derive(Debug, Clone, Copy)]
pub struct ExponentialAtmosphere {
    /// Reference density `ρ₀` at altitude `h0`.
    pub rho0: KilogramsPerCubicMeter,
    /// Reference altitude `h₀` (km).
    pub h0: Kilometers,
    /// Atmospheric scale height `H` (km).
    pub scale_height: Kilometers,
}

impl ExponentialAtmosphere {
    /// USSA-like single-layer profile representative of ~500 km LEO.
    ///
    /// Values from Vallado Table 8-4 (500 km row).  These numbers are
    /// indicative; production runs should use calibrated tables or a full
    /// atmosphere model.
    pub const LEO_500KM: Self = Self {
        rho0: KilogramsPerCubicMeter::new(6.967e-13),
        h0: Kilometers::new(500.0),
        scale_height: Kilometers::new(63.822),
    };
}

impl DensityProvider for ExponentialAtmosphere {
    fn name(&self) -> &'static str {
        "ExponentialAtmosphere"
    }

    #[inline]
    fn density(&self, altitude: Kilometers) -> Result<KilogramsPerCubicMeter, DynamicsError> {
        let h = altitude.value();
        if h < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface { altitude_km: h });
        }
        let exponent = -(h - self.h0.value()) / self.scale_height.value();
        Ok(KilogramsPerCubicMeter::new(
            self.rho0.value() * exponent.exp(),
        ))
    }
}

// =============================================================================
// ConstantDensity
// =============================================================================

/// Constant-density atmosphere.  Returns the same `ρ` for every altitude.
///
/// Useful only for controlled unit tests where the force-model logic should be
/// checked independently of the altitude–density relationship.
#[derive(Debug, Clone, Copy)]
pub struct ConstantDensity {
    /// Density returned for every altitude.
    pub rho: KilogramsPerCubicMeter,
}

impl DensityProvider for ConstantDensity {
    fn name(&self) -> &'static str {
        "ConstantDensity"
    }

    #[inline]
    fn density(&self, _altitude: Kilometers) -> Result<KilogramsPerCubicMeter, DynamicsError> {
        Ok(self.rho)
    }
}

// =============================================================================
// Nrlmsise00LiteApprox
// =============================================================================

/// Log-linear table approximation of the NRLMSISE-00 atmosphere model.
///
/// This is a **teaching / regression-test quality** density model.  It is
/// **not** suitable for operational precision-orbit determination.
///
/// ## Model description
///
/// The table below contains hand-authored altitude–density pairs for
/// F10.7 = 140 solar flux units (mean solar activity) and Ap = 15 (moderate
/// geomagnetic activity).  Density is interpolated log-linearly between table
/// entries.  Outside the 50–1000 km range the last/first bracket's slope is
/// extrapolated log-linearly.
///
/// ## Accuracy
///
/// ±50 % under mean F10.7 conditions.  No solar-flux or geomagnetic index
/// dependence is modelled.  Errors scale to a factor of 2–3 during solar
/// maximum or geomagnetic storms.
///
/// ## References
///
/// * Picone, J.M. et al. (2002), "NRLMSISE-00 empirical model of the
///   atmosphere", *J. Geophys. Res.*, 107(A12), 1468.
/// * Vallado, *Fundamentals of Astrodynamics and Applications*, Table 8-4.
#[derive(Debug, Clone, Copy, Default)]
pub struct Nrlmsise00LiteApprox;

impl DensityProvider for Nrlmsise00LiteApprox {
    fn name(&self) -> &'static str {
        "Nrlmsise00LiteApprox"
    }

    fn density(&self, altitude: Kilometers) -> Result<KilogramsPerCubicMeter, DynamicsError> {
        let h = altitude.value();
        if h < 0.0 {
            return Err(DynamicsError::AltitudeBelowSurface { altitude_km: h });
        }
        Ok(KilogramsPerCubicMeter::new(nrlmsise_interpolate(h)))
    }
}

/// Log-linear interpolation / extrapolation on [`NRLMSISE_TABLE`].
///
/// Below the table minimum (50 km) the first-point density is returned.
/// Above the table maximum (1000 km) the last bracket's log-linear slope is
/// extrapolated.
fn nrlmsise_interpolate(h_km: f64) -> f64 {
    let table = NRLMSISE_TABLE;
    let n = table.len();

    if h_km <= table[0].0 {
        return table[0].1;
    }
    if h_km >= table[n - 1].0 {
        let (h1, rho1) = table[n - 2];
        let (h2, rho2) = table[n - 1];
        let frac = (h_km - h1) / (h2 - h1);
        return (rho1.ln() + frac * (rho2.ln() - rho1.ln())).exp();
    }
    for i in 0..n - 1 {
        let (h1, rho1) = table[i];
        let (h2, rho2) = table[i + 1];
        if h_km >= h1 && h_km <= h2 {
            let frac = (h_km - h1) / (h2 - h1);
            return (rho1.ln() + frac * (rho2.ln() - rho1.ln())).exp();
        }
    }
    table[n - 1].1 // unreachable guard
}

// =============================================================================
// geodetic_altitude
// =============================================================================

/// Compute geodetic altitude above a reference ellipsoid.
///
/// Uses the iterative Bowring (1985) method; three iterations give
/// sub-millimetre convergence for altitudes up to ~36 000 km.
///
/// # Arguments
///
/// * `r` — Position in km.  **Treated as ITRF/ECEF** with the current
///   **placeholder identity GCRS↔ITRF rotation**.  Error vs. true
///   geodetic altitude is ≤ 20 km at the poles (tolerable for all
///   built-in density models).
/// * `a_eq` — Equatorial radius (km).  For WGS-84: 6 378.137 km.
/// * `f_flat` — Flattening (dimensionless).  For WGS-84: 1/298.257 223 563.
///
/// # Returns
///
/// Geodetic altitude `h` (km) above the specified ellipsoid.
///
/// # References
///
/// * Bowring, B.R. (1985), "The geodetic line and the geodetic coordinates",
///   *Survey Review*, 28(218), 276–281.
/// * Vallado (2013), *Fundamentals of Astrodynamics and Applications*, §3.4.
pub fn geodetic_altitude(
    r: &Position<Geocentric, GCRS, Kilometer>,
    a_eq: Kilometers,
    f_flat: Ratios,
) -> Kilometers {
    let x = r.x().value();
    let y = r.y().value();
    let z = r.z().value();

    let a = a_eq.value();
    let f = f_flat.value();
    let b = a * (1.0 - f);
    let e2 = 2.0 * f - f * f; // first eccentricity²
    let ep2 = (a * a - b * b) / (b * b); // second eccentricity²

    let p = (x * x + y * y).sqrt(); // equatorial distance

    // Initial parametric-latitude estimate (Bowring)
    let mut beta = (z * a).atan2(p * b);

    // 3 Bowring iterations → sub-mm convergence
    for _ in 0..3 {
        let sb3 = beta.sin().powi(3);
        let cb3 = beta.cos().powi(3);
        let phi = (z + ep2 * b * sb3).atan2(p - e2 * a * cb3);
        beta = ((1.0 - f) * phi.sin()).atan2(phi.cos());
    }

    // Final geodetic latitude
    let sb3 = beta.sin().powi(3);
    let cb3 = beta.cos().powi(3);
    let phi = (z + ep2 * b * sb3).atan2(p - e2 * a * cb3);

    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    // Normal radius of curvature
    let n_roc = a / (1.0 - e2 * sin_phi * sin_phi).sqrt();

    let h = if cos_phi.abs() > 1e-10 {
        p / cos_phi - n_roc
    } else {
        z.abs() / sin_phi.abs() - n_roc * (1.0 - e2)
    };

    Kilometers::new(h)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use affn::cartesian::Position as AffnPos;

    // ── ExponentialAtmosphere ─────────────────────────────────────────────────

    #[test]
    fn exponential_decreases_with_altitude() {
        let atm = ExponentialAtmosphere::LEO_500KM;
        assert!(
            atm.density(Kilometers::new(400.0)).unwrap().value()
                > atm.density(Kilometers::new(500.0)).unwrap().value()
        );
        assert!(
            atm.density(Kilometers::new(500.0)).unwrap().value()
                > atm.density(Kilometers::new(600.0)).unwrap().value()
        );
    }

    #[test]
    fn exponential_400km_order_of_magnitude() {
        let atm = ExponentialAtmosphere::LEO_500KM;
        let rho = atm.density(Kilometers::new(400.0)).unwrap().value();
        // LEO 400 km: ~1e-13 to 1e-10 kg/m³
        assert!(rho > 1e-13 && rho < 1e-10, "rho = {rho:.3e}");
    }

    #[test]
    fn exponential_below_surface_returns_error() {
        let atm = ExponentialAtmosphere::LEO_500KM;
        let result = atm.density(Kilometers::new(-1.0));
        assert!(
            matches!(result, Err(DynamicsError::AltitudeBelowSurface { .. })),
            "expected AltitudeBelowSurface, got {result:?}"
        );
    }

    #[test]
    fn exponential_density_is_typed() {
        let atm = ExponentialAtmosphere::LEO_500KM;
        let _typed: KilogramsPerCubicMeter = atm.density(Kilometers::new(400.0)).unwrap();
    }

    // ── ConstantDensity ───────────────────────────────────────────────────────

    #[test]
    fn constant_density_independent_of_altitude() {
        let rho = KilogramsPerCubicMeter::new(1.23e-12);
        let atm = ConstantDensity { rho };
        assert_eq!(atm.density(Kilometers::new(0.0)).unwrap().value(), 1.23e-12);
        assert_eq!(
            atm.density(Kilometers::new(1_000.0)).unwrap().value(),
            1.23e-12
        );
    }

    // ── Nrlmsise00LiteApprox ──────────────────────────────────────────────────

    #[test]
    fn nrlmsise_400km_order_of_magnitude() {
        let atm = Nrlmsise00LiteApprox;
        let rho = atm.density(Kilometers::new(400.0)).unwrap().value();
        // Table value at 400 km: 2.803e-12 kg/m³ (±50 %: ~1.4e-12 to ~4.2e-12)
        assert!(rho > 1e-12 && rho < 1e-11, "400 km rho = {rho:.3e}");
    }

    #[test]
    fn nrlmsise_decreases_with_altitude() {
        let atm = Nrlmsise00LiteApprox;
        let rho300 = atm.density(Kilometers::new(300.0)).unwrap().value();
        let rho500 = atm.density(Kilometers::new(500.0)).unwrap().value();
        let rho800 = atm.density(Kilometers::new(800.0)).unwrap().value();
        assert!(rho300 > rho500 && rho500 > rho800);
    }

    #[test]
    fn nrlmsise_below_surface_returns_error() {
        let atm = Nrlmsise00LiteApprox;
        let result = atm.density(Kilometers::new(-5.0));
        assert!(
            matches!(result, Err(DynamicsError::AltitudeBelowSurface { .. })),
            "expected AltitudeBelowSurface, got {result:?}"
        );
    }

    #[test]
    fn nrlmsise_density_is_typed() {
        let atm = Nrlmsise00LiteApprox;
        let _typed: KilogramsPerCubicMeter = atm.density(Kilometers::new(400.0)).unwrap();
    }

    #[test]
    fn nrlmsise_above_table_extrapolates() {
        let atm = Nrlmsise00LiteApprox;
        // Above 1000 km, log-linear extrapolation should give a small positive density
        let rho = atm.density(Kilometers::new(1200.0)).unwrap().value();
        assert!(rho > 0.0 && rho < 1e-14, "1200 km rho = {rho:.3e}");
    }

    // ── geodetic_altitude ─────────────────────────────────────────────────────

    #[test]
    fn geodetic_altitude_equatorial_point() {
        // A point on the equator at exactly R_earth should give h ≈ 0
        let r_eq = 6_378.137_f64;
        let pos: AffnPos<Geocentric, GCRS, Kilometer> = AffnPos::new(r_eq, 0.0_f64, 0.0_f64);
        let h = geodetic_altitude(
            &pos,
            Kilometers::new(6_378.137),
            Ratios::new(1.0 / 298.257_223_563),
        );
        assert!(h.value().abs() < 0.01, "h = {:.6}", h.value());
    }

    #[test]
    fn geodetic_altitude_leo_orbit() {
        // LEO orbit at 400 km above equator
        let r = 6_378.137 + 400.0;
        let pos: AffnPos<Geocentric, GCRS, Kilometer> = AffnPos::new(r, 0.0_f64, 0.0_f64);
        let h = geodetic_altitude(
            &pos,
            Kilometers::new(6_378.137),
            Ratios::new(1.0 / 298.257_223_563),
        );
        // At equator geocentric == geodetic → expect ≈ 400 km
        assert!((h.value() - 400.0).abs() < 0.1, "h = {:.6}", h.value());
    }
}
