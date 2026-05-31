// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ESA LISA mission — orbit reader and inter-satellite range model
//!
//! This example demonstrates how to integrate a mission-specific orbit source
//! with the generic `siderust::pod` machinery using only public library APIs.
//!
//! ## What this example shows
//!
//! 1. Parse CCSDS OEM v2.0 files (one per LISA spacecraft) with
//!    [`siderust::formats::ccsds::oem::read_oem`].
//! 2. Build a typed [`LisaEphemerisProvider`] that implements
//!    [`siderust::pod::providers::EphemerisProvider`] via cubic Hermite
//!    interpolation.
//! 3. Compute a one-way light-time-corrected inter-satellite range
//!    observation using [`InterSatRangeObs`].
//!
//! ## Running the example
//!
//! ```text
//! cargo run --example 18_lisa_pod --features pod
//! ```
//!
//! The example loads orbit files from `tests/test-data/lisa/`.
//!
//! ## References
//!
//! - Consultative Committee for Space Data Systems. (2019). Orbit Data
//!   Messages, CCSDS 502.0-B-3.
//! - Martens, W., Joffre, E. (2021). Trajectory Design for the ESA LISA
//!   Mission. *J. Astronaut. Sci.*, 68, 402–443.
//!   <https://doi.org/10.1007/s40295-021-00263-2>
//! - Danzmann, K. et al. (2017). *LISA: Laser Interferometer Space Antenna.*
//!   ESA/SRE(2017)1. <https://arxiv.org/abs/1702.00786>
#![allow(clippy::print_stdout, missing_docs, unreachable_pub)]

use std::io::Read;
use std::sync::Arc;

use affn::cartesian;
use affn::centers::{AffineCenter, ReferenceCenter};
use affn::frames::EME2000;
use qtty::unit::Kilometer;
use qtty::{Day, Meter};
use tempoch::{Time, TDB};

use siderust::formats::ccsds::oem::{read_oem, OemFile};
use siderust::formats::FormatError;
use siderust::pod::observation::error::PodObservationsError;
use siderust::pod::observation::obs_trait::{CartesianState, ObsType, Observation};
use siderust::pod::observation::provider_bundle::ProviderBundle;
use siderust::pod::providers::EphemerisProvider;
use siderust::time::JulianDate;

// ── KmPerSecond unit ─────────────────────────────────────────────────────────

type KmPerSecond = qtty::Per<Kilometer, qtty::unit::Second>;

// ── HeliocentricCenter ────────────────────────────────────────────────────────

/// Heliocentric reference centre for LISA orbit positions (Sun's centre of mass).
///
/// LISA OEM files use `CENTER_NAME = SUN`. This marker type anchors the
/// type-level centre constraint on position values.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct HeliocentricCenter;

impl ReferenceCenter for HeliocentricCenter {
    type Params = ();
    fn center_name() -> &'static str {
        "Heliocentric"
    }
}

impl AffineCenter for HeliocentricCenter {}

// ── Typed position / velocity aliases ────────────────────────────────────────

pub type LisaPosition = cartesian::Position<HeliocentricCenter, EME2000, Kilometer>;
pub type LisaVelocity = cartesian::Velocity<EME2000, KmPerSecond>;

// ── LisaSpacecraftId ─────────────────────────────────────────────────────────

/// Identifier for one of the three LISA spacecraft.
///
/// The numeric suffix matches the OEM file extension: `.oem1` → SC1, etc.
/// NAIF-style integer codes (−1001, −1002, −1003) are used for
/// [`LisaEphemerisProvider`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LisaSpacecraftId {
    SC1,
    SC2,
    SC3,
}

impl LisaSpacecraftId {
    pub fn naif_id(self) -> i32 {
        match self {
            Self::SC1 => -1001,
            Self::SC2 => -1002,
            Self::SC3 => -1003,
        }
    }

    pub fn from_naif_id(id: i32) -> Option<Self> {
        match id {
            -1001 => Some(Self::SC1),
            -1002 => Some(Self::SC2),
            -1003 => Some(Self::SC3),
            _ => None,
        }
    }
}

// ── LisaOrbitPoint ────────────────────────────────────────────────────────────

/// A single tabulated state vector from a LISA OEM file.
#[derive(Debug, Clone)]
pub struct LisaOrbitPoint {
    pub epoch: Time<TDB>,
    pub position: LisaPosition,
    pub velocity: LisaVelocity,
    pub(self) epoch_j2000_s: f64,
}

// ── LisaOrbit ─────────────────────────────────────────────────────────────────

/// Tabulated orbit for one LISA spacecraft, loaded from an OEM file.
#[derive(Debug, Clone)]
pub struct LisaOrbit {
    pub spacecraft_id: LisaSpacecraftId,
    pub points: Vec<LisaOrbitPoint>,
}

// ── LisaOrbitSet ─────────────────────────────────────────────────────────────

/// All three LISA spacecraft orbits loaded together.
#[derive(Debug, Clone)]
pub struct LisaOrbitSet {
    pub sc1: LisaOrbit,
    pub sc2: LisaOrbit,
    pub sc3: LisaOrbit,
}

impl LisaOrbitSet {
    pub fn orbit(&self, sc: LisaSpacecraftId) -> &LisaOrbit {
        match sc {
            LisaSpacecraftId::SC1 => &self.sc1,
            LisaSpacecraftId::SC2 => &self.sc2,
            LisaSpacecraftId::SC3 => &self.sc3,
        }
    }
}

// ── LisaOrbitReader ───────────────────────────────────────────────────────────

/// Reader for CCSDS OEM v2.0 LISA orbit files.
pub struct LisaOrbitReader;

impl LisaOrbitReader {
    pub fn read<R: Read>(
        reader: R,
        spacecraft_id: LisaSpacecraftId,
    ) -> Result<LisaOrbit, FormatError> {
        let oem_file = read_oem(reader)?;
        Self::convert(oem_file, spacecraft_id)
    }

    pub fn from_str(s: &str, spacecraft_id: LisaSpacecraftId) -> Result<LisaOrbit, FormatError> {
        Self::read(s.as_bytes(), spacecraft_id)
    }

    fn convert(
        oem_file: OemFile,
        spacecraft_id: LisaSpacecraftId,
    ) -> Result<LisaOrbit, FormatError> {
        let mut points: Vec<LisaOrbitPoint> = Vec::new();
        for segment in oem_file.segments {
            for state in segment.states {
                let epoch_j2000_s = jd_to_j2000_seconds(state.epoch_jd);
                let epoch = jd_to_time_tdb(state.epoch_jd)?;
                let [x, y, z] = state.position_km;
                let [vx, vy, vz] = state.velocity_km_s;
                let position = LisaPosition::new(
                    qtty::Kilometer::new(x),
                    qtty::Kilometer::new(y),
                    qtty::Kilometer::new(z),
                );
                let velocity = LisaVelocity::new(
                    qtty::Quantity::<KmPerSecond>::new(vx),
                    qtty::Quantity::<KmPerSecond>::new(vy),
                    qtty::Quantity::<KmPerSecond>::new(vz),
                );
                points.push(LisaOrbitPoint {
                    epoch,
                    epoch_j2000_s,
                    position,
                    velocity,
                });
            }
        }
        if points.is_empty() {
            return Err(FormatError::Format(
                "lisa: OEM file contains no state vectors".into(),
            ));
        }
        points.sort_by(|a, b| a.epoch_j2000_s.partial_cmp(&b.epoch_j2000_s).unwrap());
        Ok(LisaOrbit {
            spacecraft_id,
            points,
        })
    }
}

// ── LisaProviderError ─────────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum LisaProviderError {
    #[error("unknown LISA body id {0} (expected -1001, -1002, or -1003)")]
    UnknownBody(i32),
    #[error("epoch {0:.3} s (TDB J2000) is outside the covered interval [{1:.3}, {2:.3}]")]
    OutOfRange(f64, f64, f64),
}

// ── LisaEphemerisProvider ─────────────────────────────────────────────────────

/// Ephemeris provider backed by three LISA tabulated orbit files.
///
/// Implements [`EphemerisProvider`] via cubic Hermite interpolation.
///
/// `body_naif_id` mapping: −1001 → SC1, −1002 → SC2, −1003 → SC3.
pub struct LisaEphemerisProvider {
    orbits: LisaOrbitSet,
}

impl LisaEphemerisProvider {
    pub fn new(orbits: LisaOrbitSet) -> Self {
        Self { orbits }
    }
}

impl EphemerisProvider for LisaEphemerisProvider {
    type State = LisaOrbitPoint;
    type Error = LisaProviderError;

    fn state(
        &self,
        body_naif_id: i32,
        epoch_seconds_tdb: f64,
    ) -> Result<LisaOrbitPoint, LisaProviderError> {
        let sc = LisaSpacecraftId::from_naif_id(body_naif_id)
            .ok_or(LisaProviderError::UnknownBody(body_naif_id))?;
        let orbit = self.orbits.orbit(sc);
        hermite_interp(orbit, epoch_seconds_tdb)
    }
}

// ── Cubic Hermite interpolation ───────────────────────────────────────────────
//
// Given bracketing states (t₀, p₀, v₀) and (t₁, p₁, v₁) and normalised
// parameter τ = (t − t₀)/(t₁ − t₀):
//
//   h₀₀(τ) =  2τ³ − 3τ² + 1
//   h₁₀(τ) =   τ³ − 2τ² + τ   (× dt)
//   h₀₁(τ) = −2τ³ + 3τ²
//   h₁₁(τ) =   τ³ −  τ²        (× dt)
//
//   p(τ) = h₀₀·p₀ + h₁₀·dt·v₀ + h₀₁·p₁ + h₁₁·dt·v₁
//   v(τ) = dp/dt (derivative / dt)

fn hermite_interp(orbit: &LisaOrbit, t: f64) -> Result<LisaOrbitPoint, LisaProviderError> {
    let pts = &orbit.points;
    if pts.is_empty() {
        return Err(LisaProviderError::OutOfRange(t, f64::NAN, f64::NAN));
    }
    let t0 = pts.first().unwrap().epoch_j2000_s;
    let t1 = pts.last().unwrap().epoch_j2000_s;
    if t < t0 || t > t1 {
        return Err(LisaProviderError::OutOfRange(t, t0, t1));
    }
    let idx = match pts.binary_search_by(|p| {
        p.epoch_j2000_s
            .partial_cmp(&t)
            .unwrap_or(std::cmp::Ordering::Less)
    }) {
        Ok(i) => return Ok(pts[i].clone()),
        Err(i) => i.saturating_sub(1).min(pts.len() - 2),
    };

    let p0 = &pts[idx];
    let p1 = &pts[idx + 1];
    let dt = p1.epoch_j2000_s - p0.epoch_j2000_s;
    let tau = (t - p0.epoch_j2000_s) / dt;
    let tau2 = tau * tau;
    let tau3 = tau2 * tau;

    let h00 = 2.0 * tau3 - 3.0 * tau2 + 1.0;
    let h10 = tau3 - 2.0 * tau2 + tau;
    let h01 = -2.0 * tau3 + 3.0 * tau2;
    let h11 = tau3 - tau2;

    let interp_pos =
        |a: f64, da: f64, b: f64, db: f64| h00 * a + h10 * dt * da + h01 * b + h11 * dt * db;
    let x = interp_pos(
        p0.position.x().value(),
        p0.velocity.x().value(),
        p1.position.x().value(),
        p1.velocity.x().value(),
    );
    let y = interp_pos(
        p0.position.y().value(),
        p0.velocity.y().value(),
        p1.position.y().value(),
        p1.velocity.y().value(),
    );
    let z = interp_pos(
        p0.position.z().value(),
        p0.velocity.z().value(),
        p1.position.z().value(),
        p1.velocity.z().value(),
    );

    let interp_vel = |a: f64, da: f64, b: f64, db: f64| {
        let dh00 = 6.0 * tau2 - 6.0 * tau;
        let dh10 = 3.0 * tau2 - 4.0 * tau + 1.0;
        let dh01 = -6.0 * tau2 + 6.0 * tau;
        let dh11 = 3.0 * tau2 - 2.0 * tau;
        (dh00 * a + dh10 * dt * da + dh01 * b + dh11 * dt * db) / dt
    };
    let vx = interp_vel(
        p0.position.x().value(),
        p0.velocity.x().value(),
        p1.position.x().value(),
        p1.velocity.x().value(),
    );
    let vy = interp_vel(
        p0.position.y().value(),
        p0.velocity.y().value(),
        p1.position.y().value(),
        p1.velocity.y().value(),
    );
    let vz = interp_vel(
        p0.position.z().value(),
        p0.velocity.z().value(),
        p1.position.z().value(),
        p1.velocity.z().value(),
    );

    let position = LisaPosition::new(
        qtty::Kilometer::new(x),
        qtty::Kilometer::new(y),
        qtty::Kilometer::new(z),
    );
    let velocity = LisaVelocity::new(
        qtty::Quantity::<KmPerSecond>::new(vx),
        qtty::Quantity::<KmPerSecond>::new(vy),
        qtty::Quantity::<KmPerSecond>::new(vz),
    );

    let epoch_j2000_s = p0.epoch_j2000_s + tau * dt;
    let epoch_jd = j2000_seconds_to_jd(epoch_j2000_s);
    let epoch = jd_to_time_tdb(epoch_jd).unwrap_or(p0.epoch);

    Ok(LisaOrbitPoint {
        epoch,
        epoch_j2000_s,
        position,
        velocity,
    })
}

// ── Time helpers ──────────────────────────────────────────────────────────────

const J2000_JD: f64 = 2_451_545.0;

fn jd_to_j2000_seconds(jd: f64) -> f64 {
    (jd - J2000_JD) * 86_400.0
}

fn j2000_seconds_to_jd(s: f64) -> f64 {
    s / 86_400.0 + J2000_JD
}

fn jd_to_time_tdb(jd: f64) -> Result<Time<TDB>, FormatError> {
    tempoch::JulianDate::<TDB>::try_new(Day::new(jd))
        .map(Into::into)
        .map_err(|e| FormatError::Format(format!("lisa: cannot convert JD {jd} to TDB Time: {e}")))
}

// ── InterSatRangeObs ──────────────────────────────────────────────────────────

const C_KM_S: f64 = 299_792.458;

/// Inter-satellite range observation between two LISA spacecraft.
///
/// Implements one Newton iteration of light-time correction:
/// `t_emit = t_recv − ρ₀/c`.
///
/// `residual` returns `measured_m − modelled_m` (O−C residual, metres).
///
/// ## Scientific scope
///
/// LISA measures distances between spacecraft pairs via laser interferometry.
/// The measurement equation for one range observable is:
///
/// ```text
/// ρ_meas = |r_A(t_recv) − r_B(t_emit)| + ε
/// ```
///
/// where `t_emit = t_recv − ρ/c`.
#[derive(Clone)]
pub struct InterSatRangeObs {
    /// Receiving spacecraft (state at `epoch`).
    pub sc_a: LisaSpacecraftId,
    /// Transmitting spacecraft (state at `epoch − ρ/c`).
    pub sc_b: LisaSpacecraftId,
    /// Observation epoch (TT Julian date, ≈ TDB for light-time purposes).
    pub epoch: JulianDate,
    /// Measured one-way range (metres).
    pub measured_m: f64,
    /// Assumed measurement standard deviation (metres; ≈ pm-level for LISA).
    pub sigma: Meter,
    /// Shared LISA ephemeris provider.
    pub provider: Arc<LisaEphemerisProvider>,
}

impl std::fmt::Debug for InterSatRangeObs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("InterSatRangeObs")
            .field("sc_a", &self.sc_a)
            .field("sc_b", &self.sc_b)
            .field("epoch_jd", &self.epoch.value())
            .field("measured_m", &self.measured_m)
            .field("sigma", &self.sigma.value())
            .finish()
    }
}

fn jd_to_j2000_s_obs(jd: JulianDate) -> f64 {
    (jd.value() - J2000_JD) * 86_400.0
}

impl Observation for InterSatRangeObs {
    type Residual = f64;

    fn residual(
        &self,
        _state: &CartesianState,
        _providers: &dyn ProviderBundle,
    ) -> Result<f64, PodObservationsError> {
        let t_recv = jd_to_j2000_s_obs(self.epoch);
        let naif_a = self.sc_a.naif_id();
        let naif_b = self.sc_b.naif_id();

        let pt_a = self
            .provider
            .state(naif_a, t_recv)
            .map_err(|_e| PodObservationsError::LightTimeNotConverged)?;
        let pt_b0 = self
            .provider
            .state(naif_b, t_recv)
            .map_err(|_| PodObservationsError::LightTimeNotConverged)?;

        let pos_a = pt_a.position;
        let pos_b0 = pt_b0.position;
        let dx0 = pos_a.x().value() - pos_b0.x().value();
        let dy0 = pos_a.y().value() - pos_b0.y().value();
        let dz0 = pos_a.z().value() - pos_b0.z().value();
        let rho0_km = (dx0 * dx0 + dy0 * dy0 + dz0 * dz0).sqrt();

        let t_emit = t_recv - rho0_km / C_KM_S;
        let pt_b1 = self
            .provider
            .state(naif_b, t_emit)
            .map_err(|_| PodObservationsError::LightTimeNotConverged)?;

        let pos_b1 = pt_b1.position;
        let dx1 = pos_a.x().value() - pos_b1.x().value();
        let dy1 = pos_a.y().value() - pos_b1.y().value();
        let dz1 = pos_a.z().value() - pos_b1.z().value();
        let rho_km = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt();
        let rho_m = rho_km * 1_000.0;

        Ok(self.measured_m - rho_m)
    }

    fn obs_type(&self) -> ObsType {
        ObsType::InterSatRange
    }

    fn epoch(&self) -> JulianDate {
        self.epoch
    }

    fn sigma(&self) -> Meter {
        self.sigma
    }
}

// ── main ──────────────────────────────────────────────────────────────────────

fn load_orbit(path: &str, sc: LisaSpacecraftId) -> LisaOrbit {
    let raw = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot open {path}: {e}"));
    LisaOrbitReader::from_str(&raw, sc).unwrap_or_else(|e| panic!("cannot parse {path}: {e}"))
}

fn main() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test-data/lisa");

    let provider = Arc::new(LisaEphemerisProvider::new(LisaOrbitSet {
        sc1: load_orbit(
            &format!("{root}/lisa_orbit_sample.oem1"),
            LisaSpacecraftId::SC1,
        ),
        sc2: load_orbit(
            &format!("{root}/lisa_orbit_sample.oem2"),
            LisaSpacecraftId::SC2,
        ),
        sc3: load_orbit(
            &format!("{root}/lisa_orbit_sample.oem3"),
            LisaSpacecraftId::SC3,
        ),
    }));

    // Query SC1 at its first tabulated epoch.
    let sc1_t0 = provider.orbits.sc1.points[0].epoch_j2000_s;
    let pt = provider.state(-1001, sc1_t0).expect("SC1 state query");
    println!(
        "SC1 position at t₀: ({:.0}, {:.0}, {:.0}) km",
        pt.position.x().value(),
        pt.position.y().value(),
        pt.position.z().value()
    );

    // Build a synthetic inter-satellite range observation.
    let epoch_jd = j2000_seconds_to_jd(sc1_t0);
    let epoch = JulianDate::new(epoch_jd);
    let modelled_range_m = {
        let pt_a = provider.state(-1001, sc1_t0).unwrap();
        let pt_b = provider.state(-1002, sc1_t0).unwrap();
        let dx = pt_a.position.x().value() - pt_b.position.x().value();
        let dy = pt_a.position.y().value() - pt_b.position.y().value();
        let dz = pt_a.position.z().value() - pt_b.position.z().value();
        (dx * dx + dy * dy + dz * dz).sqrt() * 1_000.0
    };

    let obs = InterSatRangeObs {
        sc_a: LisaSpacecraftId::SC1,
        sc_b: LisaSpacecraftId::SC2,
        epoch,
        measured_m: modelled_range_m, // zero residual by construction
        sigma: Meter::new(1e-12),
        provider,
    };
    println!("{obs:?}");

    use siderust::astro::dynamics::{Position, Velocity};
    use siderust::coordinates::frames::GCRS;
    use siderust::pod::observation::obs_trait::CartesianState;
    use siderust::pod::observation::provider_bundle::NullProviderBundle;
    let state = CartesianState::new(
        epoch.to_j2000s(),
        Position::<GCRS>::new(0.0, 0.0, 0.0),
        Velocity::<GCRS>::new(0.0, 0.0, 0.0),
    );
    let residual = obs
        .residual(&state, &NullProviderBundle)
        .expect("range residual");
    println!("O−C residual: {residual:.6} m  (should be ~0)");
    assert!(residual.abs() < 1.0, "residual too large: {residual}");

    println!("LISA example completed successfully.");
}

#[cfg(test)]
mod tests {
    use super::*;

    const OEM_SAMPLE: &str = "\
CCSDS_OEM_VERS = 2.0\n\
CREATION_DATE  = 2024-01-01T00:00:00\n\
ORIGINATOR     = TEST\n\
\n\
META_START\n\
OBJECT_NAME          = LISA-1\n\
OBJECT_ID            = -1001\n\
CENTER_NAME          = SUN\n\
REF_FRAME            = EME2000\n\
TIME_SYSTEM          = TDB\n\
START_TIME           = 2036-02-12T12:00:00.000\n\
STOP_TIME            = 2036-02-12T13:00:00.000\n\
META_STOP\n\
\n\
2036-02-12T12:00:00.000  100000000.0  50000000.0  10000000.0  10.0  20.0  5.0\n\
2036-02-12T13:00:00.000  100036000.0  50072000.0  10018000.0  10.0  20.0  5.0\n";

    #[test]
    fn spacecraft_id_naif_roundtrip() {
        for sc in [
            LisaSpacecraftId::SC1,
            LisaSpacecraftId::SC2,
            LisaSpacecraftId::SC3,
        ] {
            assert_eq!(LisaSpacecraftId::from_naif_id(sc.naif_id()), Some(sc));
        }
        assert_eq!(LisaSpacecraftId::from_naif_id(399), None);
    }

    #[test]
    fn orbit_reader_parses_two_points() {
        let orbit = LisaOrbitReader::from_str(OEM_SAMPLE, LisaSpacecraftId::SC1).unwrap();
        assert_eq!(orbit.spacecraft_id, LisaSpacecraftId::SC1);
        assert_eq!(orbit.points.len(), 2);
        assert!((orbit.points[0].position.x().value() - 100_000_000.0).abs() < 1.0);
    }

    #[test]
    fn hermite_exact_at_tabulated_point() {
        let orbit = LisaOrbitReader::from_str(OEM_SAMPLE, LisaSpacecraftId::SC1).unwrap();
        let t0 = orbit.points[0].epoch_j2000_s;
        let pt = hermite_interp(&orbit, t0).unwrap();
        assert!((pt.position.x().value() - orbit.points[0].position.x().value()).abs() < 1e-9);
    }

    #[test]
    fn provider_out_of_range_error() {
        let orbit = LisaOrbitReader::from_str(OEM_SAMPLE, LisaSpacecraftId::SC1).unwrap();
        let t_early = orbit.points[0].epoch_j2000_s - 1.0;
        let err = hermite_interp(&orbit, t_early).unwrap_err();
        assert!(matches!(err, LisaProviderError::OutOfRange(..)));
    }

    #[test]
    fn provider_unknown_body_error() {
        let orbit = LisaOrbitReader::from_str(OEM_SAMPLE, LisaSpacecraftId::SC1).unwrap();
        let set = LisaOrbitSet {
            sc1: orbit.clone(),
            sc2: orbit.clone(),
            sc3: orbit,
        };
        let provider = LisaEphemerisProvider::new(set);
        assert!(matches!(
            provider.state(999, 0.0).unwrap_err(),
            LisaProviderError::UnknownBody(999)
        ));
    }
}
