// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Tabulated Cartesian ephemerides for POD examples and adapters.
//!
//! This module turns time-tagged Cartesian state histories into typed
//! [`EphemerisProvider`](crate::pod::providers::EphemerisProvider)
//! implementations backed by cubic Hermite interpolation.

use std::collections::HashMap;
use std::io::Read;

use affn::cartesian::{Position, Velocity};
use affn::interpolation::{CubicHermiteQuantityTable, InterpolationError, QuantityHermiteNode};
use affn::{ReferenceCenter, ReferenceFrame};
use qtty::unit::{Kilometer, Second};
use qtty::{Day, KmPerSecond, Quantity};
use tempoch::{Time, TDB};

use crate::formats::ccsds::oem::read_oem;
use crate::formats::FormatError;
use crate::pod::providers::EphemerisProvider;

const J2000_JD: f64 = 2_451_545.0;

/// A typed Cartesian state sampled from a tabulated ephemeris.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TabulatedCartesianState<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    /// State epoch in TDB.
    pub epoch: Time<TDB>,
    /// Cartesian position in kilometers.
    pub position: Position<C, F, Kilometer>,
    /// Cartesian velocity in kilometers per second.
    pub velocity: Velocity<F, KmPerSecond>,
}

/// One body's tabulated Cartesian ephemeris.
pub struct TabulatedEphemeris<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    body_naif_id: i32,
    table: CubicHermiteQuantityTable<Second, Position<C, F, Kilometer>>,
}

impl<C, F> TabulatedEphemeris<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    /// Builds an ephemeris from time-sorted states.
    pub fn new(
        body_naif_id: i32,
        states: Vec<TabulatedCartesianState<C, F>>,
    ) -> Result<Self, InterpolationError> {
        let nodes = states
            .into_iter()
            .map(|state| QuantityHermiteNode {
                x: state.epoch.to::<tempoch::J2000s>().raw(),
                value: state.position,
                derivative: state.velocity,
            })
            .collect();
        Ok(Self {
            body_naif_id,
            table: CubicHermiteQuantityTable::new(nodes)?,
        })
    }

    /// NAIF id associated with this ephemeris.
    pub fn body_naif_id(&self) -> i32 {
        self.body_naif_id
    }

    /// Returns the first and last supported epochs.
    pub fn coverage(&self) -> Result<(Time<TDB>, Time<TDB>), TabulatedEphemerisError> {
        let samples = self.table.samples();
        let start = samples[0].x.value();
        let stop = samples[samples.len() - 1].x.value();
        Ok((seconds_to_time_tdb(start)?, seconds_to_time_tdb(stop)?))
    }

    /// Evaluates this ephemeris at a typed TDB epoch.
    pub fn state_at(
        &self,
        epoch: Time<TDB>,
    ) -> Result<TabulatedCartesianState<C, F>, InterpolationError> {
        let evaluated = self.table.evaluate(epoch.to::<tempoch::J2000s>().raw())?;
        Ok(TabulatedCartesianState {
            epoch,
            position: evaluated.value,
            velocity: evaluated.derivative,
        })
    }

    /// Evaluates this ephemeris at TDB seconds since J2000.
    pub fn state_seconds(
        &self,
        epoch_seconds_tdb: f64,
    ) -> Result<TabulatedCartesianState<C, F>, TabulatedEphemerisError> {
        let epoch = seconds_to_time_tdb(epoch_seconds_tdb)?;
        self.state_at(epoch).map_err(TabulatedEphemerisError::from)
    }
}

/// Ephemeris provider backed by multiple tabulated Cartesian ephemerides.
pub struct TabulatedEphemerisProvider<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    ephemerides: HashMap<i32, TabulatedEphemeris<C, F>>,
}

impl<C, F> TabulatedEphemerisProvider<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    /// Builds a provider from one or more body ephemerides.
    pub fn new(
        ephemerides: Vec<TabulatedEphemeris<C, F>>,
    ) -> Result<Self, TabulatedEphemerisError> {
        let mut by_id = HashMap::new();
        for ephemeris in ephemerides {
            let body_naif_id = ephemeris.body_naif_id();
            if by_id.insert(body_naif_id, ephemeris).is_some() {
                return Err(TabulatedEphemerisError::DuplicateBody(body_naif_id));
            }
        }
        Ok(Self { ephemerides: by_id })
    }

    /// Returns the ephemeris for a body id, if present.
    pub fn ephemeris(&self, body_naif_id: i32) -> Option<&TabulatedEphemeris<C, F>> {
        self.ephemerides.get(&body_naif_id)
    }

    /// Evaluates one body at a typed TDB epoch.
    pub fn state_at(
        &self,
        body_naif_id: i32,
        epoch: Time<TDB>,
    ) -> Result<TabulatedCartesianState<C, F>, TabulatedEphemerisError> {
        self.ephemeris(body_naif_id)
            .ok_or(TabulatedEphemerisError::UnknownBody(body_naif_id))?
            .state_at(epoch)
            .map_err(TabulatedEphemerisError::from)
    }
}

impl<C, F> EphemerisProvider for TabulatedEphemerisProvider<C, F>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    type State = TabulatedCartesianState<C, F>;
    type Error = TabulatedEphemerisError;

    fn state(&self, body_naif_id: i32, epoch_seconds_tdb: f64) -> Result<Self::State, Self::Error> {
        self.ephemeris(body_naif_id)
            .ok_or(TabulatedEphemerisError::UnknownBody(body_naif_id))?
            .state_seconds(epoch_seconds_tdb)
    }
}

/// Errors returned by tabulated ephemeris construction and lookup.
#[derive(Debug, thiserror::Error)]
pub enum TabulatedEphemerisError {
    /// The requested body is not present.
    #[error("unknown tabulated body id {0}")]
    UnknownBody(i32),
    /// A provider was constructed with the same body id twice.
    #[error("duplicate tabulated body id {0}")]
    DuplicateBody(i32),
    /// The input epoch cannot be represented as a TDB time.
    #[error("invalid TDB epoch: {0}")]
    InvalidEpoch(String),
    /// Interpolation failed.
    #[error("interpolation failed: {0}")]
    Interpolation(#[from] InterpolationError),
    /// OEM parsing or conversion failed.
    #[error("OEM conversion failed: {0}")]
    Format(#[from] FormatError),
}

/// Reads a TDB CCSDS OEM file into a typed tabulated ephemeris.
///
/// All OEM segments must declare `TIME_SYSTEM = TDB`. Position and velocity
/// fields are interpreted with the CCSDS OEM conventional units: km and km/s.
pub fn read_oem_tdb_ephemeris<C, F, R>(
    body_naif_id: i32,
    reader: R,
) -> Result<TabulatedEphemeris<C, F>, TabulatedEphemerisError>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    R: Read,
{
    let oem_file = read_oem(reader)?;
    let mut states = Vec::new();
    for segment in oem_file.segments {
        if !segment.metadata.time_system.eq_ignore_ascii_case("TDB") {
            return Err(TabulatedEphemerisError::Format(FormatError::Format(
                "read_oem_tdb_ephemeris: expected TIME_SYSTEM = TDB".into(),
            )));
        }
        for state in segment.states {
            let epoch = jd_to_time_tdb(state.epoch_jd)?;
            let [x, y, z] = state.position_km;
            let [vx, vy, vz] = state.velocity_km_s;
            states.push(TabulatedCartesianState {
                epoch,
                position: Position::<C, F, Kilometer>::new(
                    qtty::Kilometer::new(x),
                    qtty::Kilometer::new(y),
                    qtty::Kilometer::new(z),
                ),
                velocity: Velocity::<F, KmPerSecond>::new(
                    Quantity::<KmPerSecond>::new(vx),
                    Quantity::<KmPerSecond>::new(vy),
                    Quantity::<KmPerSecond>::new(vz),
                ),
            });
        }
    }
    Ok(TabulatedEphemeris::new(body_naif_id, states)?)
}

fn jd_to_time_tdb(jd: f64) -> Result<Time<TDB>, TabulatedEphemerisError> {
    tempoch::JulianDate::<TDB>::try_new(Day::new(jd))
        .map(Into::into)
        .map_err(|e| TabulatedEphemerisError::InvalidEpoch(e.to_string()))
}

fn seconds_to_time_tdb(seconds: f64) -> Result<Time<TDB>, TabulatedEphemerisError> {
    Time::<TDB>::from_raw_j2000_seconds(qtty::Second::new(seconds))
        .map_err(|e| TabulatedEphemerisError::InvalidEpoch(e.to_string()))
}

/// Converts TDB seconds since J2000 into a Julian Date.
///
/// This helper is intended for examples and tests that need to synthesize OEM
/// rows or display provider epochs.
pub fn j2000_seconds_to_tdb_jd(seconds: qtty::Second) -> f64 {
    seconds.value() / 86_400.0 + J2000_JD
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::Heliocentric;
    use crate::coordinates::frames::EME2000;

    const OEM_SAMPLE: &str = "\
CCSDS_OEM_VERS = 2.0\n\
CREATION_DATE  = 2024-01-01T00:00:00\n\
ORIGINATOR     = TEST\n\
\n\
META_START\n\
OBJECT_NAME          = TEST\n\
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

    type LisaEphemeris = TabulatedEphemeris<Heliocentric, EME2000>;
    type LisaProvider = TabulatedEphemerisProvider<Heliocentric, EME2000>;

    #[test]
    fn reads_oem_into_typed_ephemeris() {
        let ephemeris: LisaEphemeris =
            read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        let state = ephemeris.state_seconds(epoch0_seconds()).unwrap();
        assert!((state.position.x().value() - 100_000_000.0).abs() < 1e-9);
        assert!((state.velocity.y().value() - 20.0).abs() < 1e-12);
    }

    #[test]
    fn provider_reports_unknown_body() {
        let ephemeris: LisaEphemeris =
            read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        let provider = LisaProvider::new(vec![ephemeris]).unwrap();
        assert!(matches!(
            provider.state(999, epoch0_seconds()).unwrap_err(),
            TabulatedEphemerisError::UnknownBody(999)
        ));
    }

    #[test]
    fn provider_rejects_duplicate_body() {
        let a: LisaEphemeris = read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        let b: LisaEphemeris = read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        match LisaProvider::new(vec![a, b]) {
            Err(TabulatedEphemerisError::DuplicateBody(-1001)) => {}
            Err(other) => panic!("unexpected error: {other}"),
            Ok(_) => panic!("duplicate body should fail"),
        }
    }

    #[test]
    fn provider_reports_out_of_range() {
        let ephemeris: LisaEphemeris =
            read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        let err = ephemeris.state_seconds(epoch0_seconds() - 1.0).unwrap_err();
        assert!(matches!(
            err,
            TabulatedEphemerisError::Interpolation(InterpolationError::OutOfRange { .. })
        ));
    }

    fn epoch0_seconds() -> f64 {
        let ephemeris: LisaEphemeris =
            read_oem_tdb_ephemeris(-1001, OEM_SAMPLE.as_bytes()).unwrap();
        ephemeris.table.samples()[0].x.value()
    }
}
