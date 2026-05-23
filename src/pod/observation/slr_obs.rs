// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SLR (Satellite Laser Ranging) normal-point observation model.
//!
//! ## Scientific scope
//!
//! A normal point condenses a burst of individual SLR return pulses into a
//! single compressed range value.  The two-way time-of-flight (corrected to
//! range in metres) is affected by:
//!
//! | Correction | Model |
//! |-----------|-------|
//! | Troposphere | Marini-Murray (1973) for visible-wavelength lasers |
//! | Relativistic | Schwarzschild (Shapiro) logarithmic delay |
//! | Center-of-mass | Applied externally before storing `measured_m` |
//!
//! ### Two-way range equation
//! ```text
//! P₂w = 2·ρ + ΔT_trop + ΔR_rel + ε
//! ```
//! where ρ is the one-way geometric range at mid-epoch.
//!
//! ### Marini-Murray tropospheric model
//! ```text
//! ΔR_one_way = f(λ) · A / (sin E + B/(tan E + C))
//! ```
//! with `A = 0.002277·P`, `B = 0.0014`, `C = 0.045` and
//! `f(λ)` the wavelength-dependent factor.
//!
//! ### Light-time iteration
//! One Newton iteration is performed:
//! `t_emit = t_recv − ρ/c`,
//! where `ρ` is the one-way range.  Satellite state at `t_emit` must be
//! supplied externally in `sat_pos_gcrs_km`.
//!
//! ## References
//!
//! - Marini, J. W., & Murray, C. W. (1973). *Correction of Laser Range
//!   Tracking Data for Atmospheric Refraction at Elevations Above 10 Degrees.*
//!   NASA Technical Memorandum X-70555.
//! - Petit, G., & Luzum, B. (Eds.) (2010). *IERS Conventions 2010.*
//!   IERS Technical Note 36, Chapter 9.

use crate::astro::dynamics::{Position, Velocity};
use crate::coordinates::frames::GCRS;
use crate::time::JulianDate;

use qtty::Meter;

use super::error::PodObservationsError;
use super::obs_trait::{CartesianState, ObsType, Observation};
use crate::pod::observation::provider_bundle::ProviderBundle;

/// Speed of light in m/s (exact SI).
const C_M_S: f64 = 299_792_458.0;

/// GM for Earth (EGM96/IERS 2010, m³/s²).
const GM_M3_S2: f64 = 3.986_004_418e14;

// ─── Marini-Murray troposphere ───────────────────────────────────────────────

/// Wavelength-dependent Marini-Murray factor `f(λ_μm)`.
///
/// Valid for visible/NIR pulsed-laser wavelengths (0.355–1.064 μm).
fn marini_murray_f_lambda(lambda_um: f64) -> f64 {
    0.9650 + 0.0164 / (lambda_um * lambda_um) + 0.000_228 / (lambda_um.powi(4))
}

/// One-way Marini-Murray tropospheric delay in metres.
///
/// Uses standard-atmosphere approximation: P = 1013.25 hPa, T = 15 °C.
///
/// `elevation_rad` is the elevation angle of the satellite above the local
/// horizontal (radians).  Clamped to ≥ 3° to avoid divergence.
fn marini_murray_one_way_m(elevation_rad: f64, lambda_um: f64) -> f64 {
    let e = elevation_rad.max(3.0_f64.to_radians());
    let sin_e = e.sin();
    let tan_e = e.tan();
    let f_lambda = marini_murray_f_lambda(lambda_um);
    // A = 0.002277 * P (P in hPa)
    let a = 0.002_277 * 1013.25;
    let b = 0.0014;
    let c = 0.045;
    f_lambda * a / (sin_e + b / (tan_e + c))
}

// ─── Relativistic (Schwarzschild/Shapiro) correction ────────────────────────

/// Two-way Schwarzschild (Shapiro) delay in metres for SLR.
///
/// ```text
/// ΔR_rel = 2 · (2GM/c²) · ln((r_A + r_B + ρ) / (r_A + r_B − ρ))
/// ```
///
/// Both `r_rx` and `r_sat` are geocentric distances in kilometres; `rho_km`
/// is the one-way geometric range in kilometres.
fn shapiro_slr_m(r_rx_km: f64, r_sat_km: f64, rho_km: f64) -> f64 {
    let two_gm_c2 = 2.0 * GM_M3_S2 / (C_M_S * C_M_S); // metres
    let numer = r_rx_km * 1_000.0 + r_sat_km * 1_000.0 + rho_km * 1_000.0;
    let denom = r_rx_km * 1_000.0 + r_sat_km * 1_000.0 - rho_km * 1_000.0;
    if denom > 1.0 {
        2.0 * two_gm_c2 * (numer / denom).ln()
    } else {
        0.0
    }
}

// ─── SlrNormalPointObs ───────────────────────────────────────────────────────

/// SLR normal-point observation: two-way range between a ground station and
/// an orbiting satellite.
///
/// The measurement equation is:
/// ```text
/// P₂w = 2·ρ + ΔT_trop + ΔR_rel
/// ```
///
/// `modeled_value` returns `measured_m − modelled_m` (O−C, metres).
///
/// ## Center-of-mass correction
///
/// The center-of-mass correction (satellite-specific) should be applied
/// **before** storing `measured_m`.  This struct does not apply it internally.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::slr_obs::SlrNormalPointObs;
/// use siderust::pod::observation::obs_trait::{CartesianState, Observation};
/// use siderust::pod::observation::provider_bundle::NullProviderBundle;
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// let epoch = JulianDate::new(2_451_545.0);
/// let state: CartesianState = CartesianState::new(
///     epoch.to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
/// );
///
/// let obs = SlrNormalPointObs {
///     station_id: "7840".to_string(),
///     epoch,
///     measured_m: siderust::qtty::Meters::new(0.0),
///     sigma: siderust::qtty::Meters::new(0.02),
///     station_gcrs_km: Position::<GCRS>::new(6_378.0, 0.0, 0.0),
///     sat_pos_gcrs_km: Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     sat_vel_gcrs_km_s: Velocity::<GCRS>::new(0.0, 7.5, 0.0),
///     elevation_rad: std::f64::consts::FRAC_PI_2,
///     wavelength_um: 0.532,
/// };
/// let residual = obs.modeled_value(&state, &NullProviderBundle).unwrap();
/// // Two-way range between station and satellite at 7000 km geocentric
/// // minus station at 6378 km → 622 km one-way → ~1 244 000 m two-way.
/// // Residual = 0 - modelled → should be around -1.24e6 m.
/// assert!(residual < 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SlrNormalPointObs {
    /// ILRS station identifier (e.g. `"7840"` for Herstmonceux).
    pub station_id: String,
    /// Normal-point epoch (TT Julian date, midpoint of the pass segment).
    pub epoch: JulianDate,
    /// Measured two-way range (metres, CoM-corrected).
    pub measured_m: Meter,
    /// Assumed measurement standard deviation (metres; typically 5–20 mm).
    pub sigma: Meter,
    /// Station GCRS position at epoch (km).
    pub station_gcrs_km: Position<GCRS>,
    /// Satellite GCRS position at signal-emission epoch (km).
    /// Should be iterated for light-time externally; one additional
    /// internal Newton step is applied.
    pub sat_pos_gcrs_km: Position<GCRS>,
    /// Satellite GCRS velocity at signal-emission epoch (km/s).
    pub sat_vel_gcrs_km_s: Velocity<GCRS>,
    /// Elevation angle from station to satellite (radians).
    pub elevation_rad: f64,
    /// Laser wavelength in micrometres (e.g. 0.532 for Nd:YAG green).
    pub wavelength_um: f64,
}

impl Observation for SlrNormalPointObs {
    type Residual = f64;

    fn modeled_value(
        &self,
        _state: &CartesianState,
        _providers: &dyn ProviderBundle,
    ) -> Result<f64, PodObservationsError> {
        let r_st = self.station_gcrs_km;
        let r_sat = self.sat_pos_gcrs_km;
        let v_sat = self.sat_vel_gcrs_km_s;

        // One-way geometric range
        let los = r_st - r_sat;
        let rho_km = los.magnitude().value();

        // One Newton light-time iteration: refine satellite position
        let dt_s = rho_km / (C_M_S / 1_000.0); // seconds
        let sat_pos_refined = Position::<GCRS>::new(
            r_sat.x().value() - v_sat.x().value() * dt_s,
            r_sat.y().value() - v_sat.y().value() * dt_s,
            r_sat.z().value() - v_sat.z().value() * dt_s,
        );
        let los2 = r_st - sat_pos_refined;
        let rho_km = los2.magnitude().value();
        let rho_m = rho_km * 1_000.0;

        // Two-way geometric range
        let two_way_m = 2.0 * rho_m;

        // Marini-Murray troposphere (two-way)
        let trop_m = 2.0 * marini_murray_one_way_m(self.elevation_rad, self.wavelength_um);

        // Shapiro relativistic correction (two-way SLR)
        let r_rx_km = {
            let x = r_st.x().value();
            let y = r_st.y().value();
            let z = r_st.z().value();
            (x * x + y * y + z * z).sqrt()
        };
        let r_sat_km = {
            let x = sat_pos_refined.x().value();
            let y = sat_pos_refined.y().value();
            let z = sat_pos_refined.z().value();
            (x * x + y * y + z * z).sqrt()
        };
        let shapiro_m = shapiro_slr_m(r_rx_km, r_sat_km, rho_km);

        let modelled = two_way_m + trop_m + shapiro_m;
        Ok(self.measured_m.value() - modelled)
    }

    fn obs_type(&self) -> ObsType {
        ObsType::SlrNormalPoint
    }

    fn epoch(&self) -> JulianDate {
        self.epoch
    }

    fn sigma(&self) -> Meter {
        self.sigma
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::pod::observation::provider_bundle::NullProviderBundle;
    use approx::assert_abs_diff_eq;
    use qtty::Meter;

    const EPOCH: fn() -> JulianDate = || JulianDate::new(2_451_545.0);

    fn make_state(pos: Position<GCRS>) -> CartesianState {
        CartesianState::new(
            EPOCH().to_j2000s(),
            pos,
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn two_way_range_exact_for_consistent_obs() {
        let sat_pos = Position::<GCRS>::new(7_000.0, 0.0, 0.0);
        let sta_pos = Position::<GCRS>::new(6_378.0, 0.0, 0.0);
        let state = make_state(sat_pos);

        // Probe: compute modelled value, then set measured = modelled
        let obs_probe = SlrNormalPointObs {
            station_id: "7840".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(0.02),
            station_gcrs_km: sta_pos,
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: Velocity::<GCRS>::new(0.0, 7.5, 0.0),
            elevation_rad: std::f64::consts::FRAC_PI_2,
            wavelength_um: 0.532,
        };
        let neg_res = obs_probe
            .modeled_value(&state, &NullProviderBundle)
            .unwrap();
        let modelled = -neg_res;

        let obs = SlrNormalPointObs {
            measured_m: Meter::new(modelled),
            ..obs_probe
        };
        let r = obs.modeled_value(&state, &NullProviderBundle).unwrap();
        assert_abs_diff_eq!(r, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn zenith_trop_delay_near_5m_two_way() {
        let d = 2.0 * marini_murray_one_way_m(std::f64::consts::FRAC_PI_2, 0.532);
        // Two-way zenith delay ≈ 4–6 m for standard atmosphere
        assert!(d > 4.0 && d < 6.0, "two-way zenith trop = {d:.3} m");
    }

    #[test]
    fn marini_murray_wavelength_factor_532nm() {
        let f = marini_murray_f_lambda(0.532);
        // Expected ≈ 1.026
        assert!(f > 1.02 && f < 1.04, "f(532 nm) = {f:.4}");
    }

    #[test]
    fn shapiro_delay_positive_and_small() {
        // LEO satellite at ~7000 km, station at ~6378 km, range ~622 km
        let d = shapiro_slr_m(6_378.0, 7_000.0, 622.0);
        assert!(d > 0.0 && d < 0.1, "Shapiro = {d:.6} m");
    }
}
