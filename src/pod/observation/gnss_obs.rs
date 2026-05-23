// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! GNSS pseudorange and carrier-phase observation models.
//!
//! ## Scientific scope
//!
//! ### Pseudorange
//! `P = ρ + c·(dt_r − dt_s) + T + I + ε`
//!
//! | Term | Source |
//! |------|--------|
//! | ρ    | Geometric range via typed `affn` position difference |
//! | c·dt_r | Receiver clock bias from [`ProviderBundle`] |
//! | c·dt_s | Satellite clock bias from [`ProviderBundle`] |
//! | T    | Saastamoinen tropospheric model (dry + wet, elevation-dependent) |
//! | I    | Klobuchar or iono-free combination |
//! | Sagnac | Earth-rotation correction (LEO/MEO regime) |
//! | Rel  | Relativistic GPS clock term (−2·r·v/c) |
//!
//! ### Carrier phase
//! `L = ρ + c·(dt_r − dt_s) + T − I + λ·N + ε`
//!
//! Same geometric and atmospheric models as pseudorange; ionosphere enters
//! with the opposite sign.  The integer ambiguity `N` is stored as `i32`
//! and resolved externally.
//!
//! ## References
//!
//! - Misra, P., & Enge, P. (2012). *Global Positioning System: Signals,
//!   Measurements, and Performance* (2nd ed.). Ganga-Jamuna Press.
//! - Saastamoinen, J. (1972). Atmospheric correction for the troposphere and
//!   stratosphere in radio ranging satellites.
//!   *Geophysical Monograph Series*, 15, 247–251.
//! - Klobuchar, J. A. (1987). Ionospheric time-delay algorithm for
//!   single-frequency GPS users.
//!   *IEEE Transactions on Aerospace and Electronic Systems*, 23(3), 325–331.

use crate::astro::dynamics::forces::OMEGA_EARTH_RAD_S;
use crate::astro::dynamics::{Position, Velocity};
use crate::coordinates::frames::GCRS;
use crate::time::JulianDate;

use qtty::unit::{Kilometer, Second};
use qtty::velocity::C;
use qtty::{Hertz, Meter, Per};

use super::error::PodObservationsError;
use super::obs_trait::{CartesianState, ObsType, Observation, PhaseResidual};
use crate::pod::observation::provider_bundle::ProviderBundle;

// ─── Speed of light ──────────────────────────────────────────────────────────

/// Speed of light in m/s (exact SI definition).
const C_M_S: f64 = 299_792_458.0;

// ─── Ionosphere model ────────────────────────────────────────────────────────

/// Parameters for the Klobuchar single-frequency ionospheric model.
///
/// Coefficients are broadcast in the GPS navigation message.
/// Receiver position and geometry must be supplied for each epoch.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::gnss_obs::KlobucharParams;
///
/// let p = KlobucharParams::gps_default();
/// assert_eq!(p.alpha[0], 2.0e-8);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct KlobucharParams {
    /// GPS α coefficients [s, s/sc, s/sc², s/sc³].
    pub alpha: [f64; 4],
    /// GPS β coefficients [s, s/sc, s/sc², s/sc³].
    pub beta: [f64; 4],
    /// Receiver geodetic latitude in semi-circles (−0.5 to +0.5).
    pub receiver_lat_sc: f64,
    /// Receiver longitude in semi-circles (0 to 2).
    pub receiver_lon_sc: f64,
    /// Elevation angle to satellite in semi-circles (0 to 0.5).
    pub elevation_sc: f64,
    /// Azimuth angle from receiver to satellite in semi-circles (0 to 2).
    pub azimuth_sc: f64,
    /// UTC time of day in seconds.
    pub utc_seconds: f64,
}

impl KlobucharParams {
    /// GPS default coefficients (zero-deviation ionosphere).
    pub fn gps_default() -> Self {
        Self {
            alpha: [2.0e-8, 1.490e-8, -5.960e-8, -5.960e-8],
            beta: [90_112.0, 65_536.0, -131_072.0, -131_072.0],
            receiver_lat_sc: 0.0,
            receiver_lon_sc: 0.0,
            elevation_sc: 0.25,
            azimuth_sc: 0.0,
            utc_seconds: 0.0,
        }
    }
}

/// Ionospheric delay model selection.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::gnss_obs::{IonoModel, KlobucharParams};
///
/// let _iono_free = IonoModel::IonoFree;
/// let _klob = IonoModel::Klobuchar(KlobucharParams::gps_default());
/// let _none = IonoModel::None;
/// ```
#[derive(Debug, Clone)]
pub enum IonoModel {
    /// Iono-free linear combination already applied to the observable; `I = 0`.
    IonoFree,
    /// Klobuchar broadcast model (single-frequency users).
    Klobuchar(KlobucharParams),
    /// No ionospheric correction (single-frequency, no model available).
    None,
}

/// Tropospheric delay model selection.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::gnss_obs::TropModel;
///
/// let _saas = TropModel::Saastamoinen { elevation_rad: std::f64::consts::FRAC_PI_2 };
/// let _none = TropModel::None;
/// ```
#[derive(Debug, Clone, Copy)]
pub enum TropModel {
    /// Saastamoinen standard-atmosphere model.  The elevation angle (rad) is
    /// the angle from the local horizontal to the satellite.
    Saastamoinen {
        /// Elevation angle in radians (0 = horizontal, π/2 = zenith).
        elevation_rad: f64,
    },
    /// No tropospheric correction (applicable for LEO receivers above the
    /// troposphere or when a pre-corrected observable is supplied).
    None,
}

// ─── Atmosphere models ───────────────────────────────────────────────────────

/// Saastamoinen tropospheric delay (dry + wet) in metres.
///
/// Uses standard-atmosphere values: P = 1013.25 hPa, T = 288.15 K,
/// e = 11.691 hPa (≈ 50 % relative humidity at sea level).
///
/// Reference: Saastamoinen (1972).
fn saastamoinen_m(elevation_rad: f64) -> f64 {
    let sin_e = elevation_rad.sin().max(0.07); // clamp to avoid blow-up below ~4°
                                               // Dry: 0.002277 * P / sin(E)
    let dry = 0.002_277 * 1013.25 / sin_e;
    // Wet: 0.002277 * (1255/T + 0.05) * e / sin(E)
    let wet = 0.002_277 * (1255.0 / 288.15 + 0.05) * 11.691 / sin_e;
    dry + wet
}

/// Klobuchar L1 ionospheric delay in metres.
///
/// Reference: Klobuchar (1987), IS-GPS-200 section 20.3.3.5.
fn klobuchar_m(p: &KlobucharParams) -> f64 {
    let e = p.elevation_sc; // elevation in semi-circles

    // Earth central angle (semi-circles)
    let psi = 0.0137 / (e + 0.11) - 0.022;

    // Sub-ionospheric latitude (semi-circles), clamped to ±0.416 sc
    let phi_i = (p.receiver_lat_sc + psi * (p.azimuth_sc * std::f64::consts::PI).cos())
        .clamp(-0.416, 0.416);

    // Sub-ionospheric longitude (semi-circles)
    let cos_phi = (phi_i * std::f64::consts::PI).cos();
    let lambda_i = if cos_phi.abs() > 1e-9 {
        p.receiver_lon_sc + psi * (p.azimuth_sc * std::f64::consts::PI).sin() / cos_phi
    } else {
        p.receiver_lon_sc
    };

    // Geomagnetic latitude (semi-circles)
    let phi_m = phi_i + 0.064 * ((lambda_i - 1.617) * std::f64::consts::PI).cos();

    // Local time (seconds)
    let mut t = 4.32e4 * lambda_i + p.utc_seconds;
    t -= (t / 86400.0).floor() * 86400.0; // wrap to [0, 86400)

    // Period (seconds)
    let period_raw: f64 = p.beta[0]
        + p.beta[1] * phi_m
        + p.beta[2] * phi_m * phi_m
        + p.beta[3] * phi_m * phi_m * phi_m;
    let period = period_raw.max(72_000.0);

    // Amplitude (seconds)
    let amp_raw: f64 = p.alpha[0]
        + p.alpha[1] * phi_m
        + p.alpha[2] * phi_m * phi_m
        + p.alpha[3] * phi_m * phi_m * phi_m;
    let amp = amp_raw.max(0.0);

    // Obliquity factor
    let obliquity = 1.0 + 16.0 * (0.53 - e).powi(3);

    // Phase argument
    let x = 2.0 * std::f64::consts::PI * (t - 50_400.0) / period;

    // Vertical delay (seconds)
    let i_v = if x.abs() < std::f64::consts::FRAC_PI_2 {
        5.0e-9 + amp * x.cos()
    } else {
        5.0e-9
    };

    // Slant delay (metres)
    obliquity * i_v * C_M_S
}

// ─── Internal range helpers ──────────────────────────────────────────────────

/// Sagnac (Earth-rotation) range correction in km.
fn sagnac_km(gps_pos_km: Position<GCRS>, rx_pos_km: Position<GCRS>) -> f64 {
    let c_km_s = C.to::<Per<Kilometer, Second>>().value();
    OMEGA_EARTH_RAD_S.value()
        * (gps_pos_km.x().value() * rx_pos_km.y().value()
            - gps_pos_km.y().value() * rx_pos_km.x().value())
        / c_km_s
}

/// Relativistic GPS satellite clock correction in metres.
fn relativistic_gps_clock_m(pos_km: Position<GCRS>, vel_km_s: Velocity<GCRS>) -> f64 {
    let dot = pos_km.x().value() * vel_km_s.x().value()
        + pos_km.y().value() * vel_km_s.y().value()
        + pos_km.z().value() * vel_km_s.z().value();
    -2.0 * dot * 1_000.0 / C.value()
}

/// Compute geometric range (m) and receiver-to-satellite unit vector.
fn geometric_range_m(state: &CartesianState, sat_pos_km: Position<GCRS>) -> (f64, [f64; 3]) {
    let los = state.position - sat_pos_km;
    let rho_km = los.magnitude().value();
    let rho_m = rho_km * 1_000.0;
    let hat = if rho_km > 0.0 {
        [
            los.x().value() / rho_km,
            los.y().value() / rho_km,
            los.z().value() / rho_km,
        ]
    } else {
        [0.0; 3]
    };
    (rho_m, hat)
}

// ─── GnssPseudorangeObs ──────────────────────────────────────────────────────

/// GNSS pseudorange observation between a spacecraft receiver and one GNSS
/// transmitter.
///
/// The measurement equation is:
/// ```text
/// P = ρ + c·(dt_r − dt_s) + T + I + Sagnac + Rel + ε
/// ```
///
/// `modeled_value` returns `measured_m − modelled_m` (O−C residual, metres).
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::gnss_obs::{GnssPseudorangeObs, IonoModel, TropModel};
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
/// let sat_pos = Position::<GCRS>::new(26_560.0, 0.0, 0.0);
/// let sat_vel = Velocity::<GCRS>::new(0.0, 3.87, 0.0);
///
/// // Build a pseudorange obs with known geometry
/// let obs = GnssPseudorangeObs {
///     prn: "G01".to_string(),
///     epoch,
///     measured_m: siderust::qtty::Meters::new(19_560_000.0),
///     sigma: siderust::qtty::Meters::new(1.0),
///     sat_pos_gcrs_km: sat_pos,
///     sat_vel_gcrs_km_s: sat_vel,
///     trop: TropModel::None,
///     iono: IonoModel::IonoFree,
///     frequency_hz: siderust::qtty::Hertzs::new(1_575_420_000.0),
/// };
/// let residual = obs.modeled_value(&state, &NullProviderBundle).unwrap();
/// // Residual should be small (depends on exact geometry + Sagnac + Rel)
/// assert!(residual.abs() < 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct GnssPseudorangeObs {
    /// GNSS satellite PRN identifier (e.g. `"G01"`).
    pub prn: String,
    /// Observation epoch (TT Julian date).
    pub epoch: JulianDate,
    /// Measured pseudorange (metres).
    pub measured_m: Meter,
    /// Assumed measurement standard deviation (metres).
    pub sigma: Meter,
    /// GNSS satellite GCRS position at signal-emission epoch (km).
    pub sat_pos_gcrs_km: Position<GCRS>,
    /// GNSS satellite GCRS velocity at signal-emission epoch (km/s).
    pub sat_vel_gcrs_km_s: Velocity<GCRS>,
    /// Tropospheric delay model.
    pub trop: TropModel,
    /// Ionospheric delay model.
    pub iono: IonoModel,
    /// Signal carrier frequency (used when Klobuchar I ≠ zero).
    pub frequency_hz: Hertz,
}

impl Observation for GnssPseudorangeObs {
    type Residual = f64;

    fn modeled_value(
        &self,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Result<f64, PodObservationsError> {
        let sat_pos = self.sat_pos_gcrs_km;
        let sat_vel = self.sat_vel_gcrs_km_s;

        let (rho_m, _) = geometric_range_m(state, sat_pos);
        let sagnac_m = sagnac_km(sat_pos, state.position) * 1_000.0;
        let rel_m = relativistic_gps_clock_m(sat_pos, sat_vel);

        let trop_m = match self.trop {
            TropModel::None => 0.0,
            TropModel::Saastamoinen { elevation_rad } => saastamoinen_m(elevation_rad),
        };

        let f_hz = self.frequency_hz.value();
        let iono_m = match &self.iono {
            IonoModel::IonoFree => 0.0,
            IonoModel::None => 0.0,
            IonoModel::Klobuchar(p) => {
                let iono_l1 = klobuchar_m(p);
                // Scale to other frequencies: I ∝ 1/f²
                let f_l1_hz = 1_575_420_000.0_f64;
                if f_hz > 0.0 && (f_hz - f_l1_hz).abs() > 1e6 {
                    iono_l1 * (f_l1_hz / f_hz).powi(2)
                } else {
                    iono_l1
                }
            }
        };

        let dt_r = providers.receiver_clock_m().value();
        let dt_s = providers
            .gnss_satellite_clock_m(&self.prn, self.epoch)
            .value();

        let modelled = rho_m + dt_r - dt_s + trop_m + iono_m + sagnac_m + rel_m;
        Ok(self.measured_m.value() - modelled)
    }

    fn obs_type(&self) -> ObsType {
        ObsType::GnssPseudorange
    }

    fn epoch(&self) -> JulianDate {
        self.epoch
    }

    fn sigma(&self) -> Meter {
        self.sigma
    }
}

// ─── GnssCarrierPhaseObs ─────────────────────────────────────────────────────

/// GNSS carrier-phase observation between a spacecraft receiver and one GNSS
/// transmitter.
///
/// The measurement equation is:
/// ```text
/// L = ρ + c·(dt_r − dt_s) + T − I + λ·N + ε
/// ```
///
/// `modeled_value` returns a [`PhaseResidual`] carrying the O−C in metres and
/// in carrier cycles.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::gnss_obs::{GnssCarrierPhaseObs, IonoModel, TropModel};
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
/// let obs = GnssCarrierPhaseObs {
///     prn: "G01".to_string(),
///     epoch,
///     measured_m: siderust::qtty::Meters::new(19_560_000.0),
///     sigma: siderust::qtty::Meters::new(0.003),
///     sat_pos_gcrs_km: Position::<GCRS>::new(26_560.0, 0.0, 0.0),
///     sat_vel_gcrs_km_s: Velocity::<GCRS>::new(0.0, 3.87, 0.0),
///     trop: TropModel::None,
///     iono: IonoModel::IonoFree,
///     frequency_hz: siderust::qtty::Hertzs::new(1_575_420_000.0),
///     integer_ambiguity: 0,
/// };
/// let r = obs.modeled_value(&state, &NullProviderBundle).unwrap();
/// assert!(r.residual_m.abs() < 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct GnssCarrierPhaseObs {
    /// GNSS satellite PRN identifier.
    pub prn: String,
    /// Observation epoch (TT Julian date).
    pub epoch: JulianDate,
    /// Measured carrier-phase range (metres, scaled by wavelength externally).
    pub measured_m: Meter,
    /// Assumed measurement standard deviation (metres; typically ≈ λ/100).
    pub sigma: Meter,
    /// GNSS satellite GCRS position at signal-emission epoch (km).
    pub sat_pos_gcrs_km: Position<GCRS>,
    /// GNSS satellite GCRS velocity at signal-emission epoch (km/s).
    pub sat_vel_gcrs_km_s: Velocity<GCRS>,
    /// Tropospheric delay model.
    pub trop: TropModel,
    /// Ionospheric delay model (enters with opposite sign vs pseudorange).
    pub iono: IonoModel,
    /// Signal carrier frequency (used to compute wavelength λ = c/f).
    pub frequency_hz: Hertz,
    /// Integer phase ambiguity N (resolved externally; set to 0 for float).
    pub integer_ambiguity: i32,
}

impl Observation for GnssCarrierPhaseObs {
    type Residual = PhaseResidual;

    fn modeled_value(
        &self,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Result<PhaseResidual, PodObservationsError> {
        let sat_pos = self.sat_pos_gcrs_km;
        let sat_vel = self.sat_vel_gcrs_km_s;

        let (rho_m, _) = geometric_range_m(state, sat_pos);
        let sagnac_m = sagnac_km(sat_pos, state.position) * 1_000.0;
        let rel_m = relativistic_gps_clock_m(sat_pos, sat_vel);

        let trop_m = match self.trop {
            TropModel::None => 0.0,
            TropModel::Saastamoinen { elevation_rad } => saastamoinen_m(elevation_rad),
        };

        let f_hz = self.frequency_hz.value();
        let iono_m = match &self.iono {
            IonoModel::IonoFree => 0.0,
            IonoModel::None => 0.0,
            IonoModel::Klobuchar(p) => {
                let iono_l1 = klobuchar_m(p);
                let f_l1_hz = 1_575_420_000.0_f64;
                if f_hz > 0.0 && (f_hz - f_l1_hz).abs() > 1e6 {
                    iono_l1 * (f_l1_hz / f_hz).powi(2)
                } else {
                    iono_l1
                }
            }
        };

        let dt_r = providers.receiver_clock_m().value();
        let dt_s = providers
            .gnss_satellite_clock_m(&self.prn, self.epoch)
            .value();

        let wavelength_m = if f_hz > 0.0 {
            C_M_S / f_hz
        } else {
            0.1903 // GPS L1 ≈ 19.03 cm
        };
        let ambiguity_m = self.integer_ambiguity as f64 * wavelength_m;

        // Phase: L = ρ + c·(dt_r - dt_s) + T - I + λN
        let modelled = rho_m + dt_r - dt_s + trop_m - iono_m + sagnac_m + rel_m + ambiguity_m;
        let residual_m = self.measured_m.value() - modelled;

        Ok(PhaseResidual {
            residual_m,
            cycles: residual_m / wavelength_m,
        })
    }

    fn obs_type(&self) -> ObsType {
        ObsType::GnssCarrierPhase
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
    use qtty::{Hertz, Meter};

    const EPOCH: fn() -> JulianDate = || JulianDate::new(2_451_545.0);

    fn leo_state() -> CartesianState {
        CartesianState::new(
            EPOCH().to_j2000s(),
            Position::<GCRS>::new(7_000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn gps_sat() -> (Position<GCRS>, Velocity<GCRS>) {
        (
            Position::<GCRS>::new(26_560.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 3.87, 0.0),
        )
    }

    #[test]
    fn pseudorange_residual_near_zero_for_consistent_obs() {
        let (sat_pos, sat_vel) = gps_sat();
        let state = leo_state();
        // Build an obs whose measured_m equals the computed modelled value.
        let obs_probe = GnssPseudorangeObs {
            prn: "G01".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(1.0),
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: sat_vel,
            trop: TropModel::None,
            iono: IonoModel::IonoFree,
            frequency_hz: Hertz::new(1_575_420_000.0),
        };
        // First get the modelled value
        let neg_residual = obs_probe
            .modeled_value(&state, &NullProviderBundle)
            .unwrap();
        // neg_residual = 0 - modelled  → modelled = -neg_residual
        let modelled = -neg_residual;

        let obs = GnssPseudorangeObs {
            measured_m: Meter::new(modelled),
            ..obs_probe
        };
        let r = obs.modeled_value(&state, &NullProviderBundle).unwrap();
        assert_abs_diff_eq!(r, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn pseudorange_positive_trop_increases_range() {
        let (sat_pos, sat_vel) = gps_sat();
        let state = leo_state();
        let make = |trop| GnssPseudorangeObs {
            prn: "G01".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(1.0),
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: sat_vel,
            trop,
            iono: IonoModel::IonoFree,
            frequency_hz: Hertz::new(1_575_420_000.0),
        };
        let r_none = make(TropModel::None)
            .modeled_value(&state, &NullProviderBundle)
            .unwrap();
        let r_trop = make(TropModel::Saastamoinen {
            elevation_rad: std::f64::consts::FRAC_PI_2,
        })
        .modeled_value(&state, &NullProviderBundle)
        .unwrap();
        // Troposphere adds to the range → modelled_trop > modelled_none →
        // residual_trop (= 0 − modelled) is more negative than residual_none.
        assert!(r_trop < r_none, "trop should increase modelled range");
    }

    #[test]
    fn klobuchar_delay_positive() {
        let delay = klobuchar_m(&KlobucharParams::gps_default());
        assert!(delay > 0.0 && delay < 100.0, "klobuchar delay={delay:.3}m");
    }

    #[test]
    fn saastamoinen_zenith_approx_2p4m() {
        let d = saastamoinen_m(std::f64::consts::FRAC_PI_2);
        // Zenith delay ~2.4 m for standard atmosphere
        assert!(d > 2.0 && d < 3.0, "zenith trop={d:.3}m");
    }

    #[test]
    fn carrier_phase_iono_opposite_sign_to_pseudorange() {
        let (sat_pos, sat_vel) = gps_sat();
        let state = leo_state();
        let klob = IonoModel::Klobuchar(KlobucharParams::gps_default());

        let pr = GnssPseudorangeObs {
            prn: "G01".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(1.0),
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: sat_vel,
            trop: TropModel::None,
            iono: klob.clone(),
            frequency_hz: Hertz::new(1_575_420_000.0),
        };
        let cp = GnssCarrierPhaseObs {
            prn: "G01".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(0.003),
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: sat_vel,
            trop: TropModel::None,
            iono: klob,
            frequency_hz: Hertz::new(1_575_420_000.0),
            integer_ambiguity: 0,
        };

        let r_pr = pr.modeled_value(&state, &NullProviderBundle).unwrap();
        let r_cp = cp.modeled_value(&state, &NullProviderBundle).unwrap();

        // Pseudorange: residual = 0 - (rho + I + ...) → negative
        // Carrier:     residual = 0 - (rho - I + ...) → less negative than pseudorange
        // Thus r_cp > r_pr (less negative)
        assert!(
            r_cp.residual_m > r_pr,
            "carrier more negative than expected"
        );
    }

    #[test]
    fn carrier_phase_residual_in_cycles_consistent() {
        let (sat_pos, sat_vel) = gps_sat();
        let state = leo_state();
        let obs = GnssCarrierPhaseObs {
            prn: "G01".to_string(),
            epoch: EPOCH(),
            measured_m: Meter::new(0.0),
            sigma: Meter::new(0.003),
            sat_pos_gcrs_km: sat_pos,
            sat_vel_gcrs_km_s: sat_vel,
            trop: TropModel::None,
            iono: IonoModel::IonoFree,
            frequency_hz: Hertz::new(1_575_420_000.0),
            integer_ambiguity: 0,
        };
        let r = obs.modeled_value(&state, &NullProviderBundle).unwrap();
        let wavelength = C_M_S / 1_575_420_000.0;
        assert_abs_diff_eq!(r.cycles, r.residual_m / wavelength, epsilon = 1e-9);
    }
}
