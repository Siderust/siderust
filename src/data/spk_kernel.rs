// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Valles Puig, Ramon

//! Runtime SPK kernel stacks for body-center and barycenter states.

use super::spk::{self, IndexedSegmentData};
use super::DataError;
use crate::calculus::ephemeris::{EphemerisError, MajorPlanet, PlanetPoint};
use crate::calculus::jpl::eval::DynSegmentDescriptor;
use crate::coordinates::frames::ICRF;
use crate::qtty::{Kilometer, Kilometers};
use crate::time::{JulianDate, TDB};
use affn::Displacement;
use std::path::Path;

const SSB_ID: i32 = 0;
const EMB_ID: i32 = 3;
const EARTH_ID: i32 = 399;
const MOON_ID: i32 = 301;
const EARTH_MOON_RATIO: f64 = 81.300_569_074_190_62;
const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;

type IcrfKm = Displacement<ICRF, Kilometer>;

/// Error returned while resolving a runtime SPK target-center chain.
#[derive(Debug)]
pub enum SpkKernelError {
    /// No loaded segment can resolve the requested NAIF target.
    MissingSegment {
        /// NAIF target ID that could not be resolved.
        target_id: i32,
    },
    /// Loaded segments form a target-center cycle.
    CenterCycle {
        /// NAIF target ID at which the cycle was detected.
        target_id: i32,
    },
    /// A Chebyshev segment does not cover the requested epoch.
    Ephemeris(EphemerisError),
}

impl core::fmt::Display for SpkKernelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingSegment { target_id } => {
                write!(f, "no loaded SPK segment resolves NAIF target {target_id}")
            }
            Self::CenterCycle { target_id } => {
                write!(
                    f,
                    "SPK target-center chain contains a cycle at NAIF target {target_id}"
                )
            }
            Self::Ephemeris(err) => err.fmt(f),
        }
    }
}

impl std::error::Error for SpkKernelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Ephemeris(err) => Some(err),
            _ => None,
        }
    }
}

impl From<EphemerisError> for SpkKernelError {
    fn from(value: EphemerisError) -> Self {
        Self::Ephemeris(value)
    }
}

struct KernelSegment {
    target_id: i32,
    center_id: i32,
    start_et: f64,
    end_et: f64,
    descriptor: DynSegmentDescriptor,
}

impl KernelSegment {
    fn from_indexed(segment: IndexedSegmentData) -> Self {
        Self {
            target_id: segment.target_id,
            center_id: segment.center_id,
            start_et: segment.start_et,
            end_et: segment.end_et,
            descriptor: DynSegmentDescriptor::from_spk(&segment.data),
        }
    }

    fn covers(&self, jd: JulianDate) -> bool {
        let jd_tdb = jd.to_scale::<TDB>().raw().value();
        let et = (jd_tdb - 2_451_545.0) * SECONDS_PER_DAY;
        et.is_finite() && et >= self.start_et && et <= self.end_et
    }

    fn out_of_range(&self, jd: JulianDate) -> SpkKernelError {
        SpkKernelError::Ephemeris(EphemerisError::OutOfRange {
            jd: jd.to_scale::<TDB>().raw().value(),
            start_jd: 2_451_545.0 + self.start_et / SECONDS_PER_DAY,
            end_jd: 2_451_545.0 + self.end_et / SECONDS_PER_DAY,
        })
    }

    fn try_position(&self, jd: JulianDate) -> Result<IcrfKm, EphemerisError> {
        self.descriptor.try_position(jd.to_scale::<TDB>())
    }
}

/// Runtime-loaded stack of J2000 SPK Type 2/3 segments.
///
/// Later paths take precedence when more than one loaded segment covers the
/// same target. Each state request resolves a target-center chain back to the
/// solar-system barycenter before subtracting the requested observer.
#[derive(Default)]
pub struct SpkKernelSet {
    segments: Vec<KernelSegment>,
}

impl SpkKernelSet {
    /// Parse and stack one or more BSP/SPK files.
    pub fn from_paths<I, P>(paths: I) -> Result<Self, DataError>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        let mut set = Self::default();
        for path in paths {
            let bytes = std::fs::read(path.as_ref())?;
            set.push_bytes(&bytes)?;
        }
        if set.segments.is_empty() {
            return Err(DataError::Parse(
                "SPK kernel set contains no supported segments".to_string(),
            ));
        }
        Ok(set)
    }

    /// Add all supported segments from an in-memory BSP/SPK.
    pub fn push_bytes(&mut self, bytes: &[u8]) -> Result<(), DataError> {
        self.segments.extend(
            spk::parse_indexed_segments(bytes)?
                .into_iter()
                .map(KernelSegment::from_indexed),
        );
        Ok(())
    }

    /// Number of indexed SPK segments in this kernel stack.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Resolve a geometric target-minus-observer state in ICRF kilometres.
    pub fn try_geometric_state(
        &self,
        target_naif_id: i32,
        observer_naif_id: i32,
        jd: JulianDate,
    ) -> Result<IcrfKm, SpkKernelError> {
        let target = self.try_position_from_ssb(target_naif_id, jd, &mut Vec::new())?;
        let observer = self.try_position_from_ssb(observer_naif_id, jd, &mut Vec::new())?;
        Ok(target - observer)
    }

    /// Resolve a major-planet point relative to Earth's center.
    pub fn try_major_planet_geocentric(
        &self,
        planet: MajorPlanet,
        point: PlanetPoint,
        jd: JulianDate,
    ) -> Result<IcrfKm, SpkKernelError> {
        self.try_geometric_state(planet.naif_id(point), EARTH_ID, jd)
    }

    fn try_position_from_ssb(
        &self,
        target_id: i32,
        jd: JulianDate,
        visited: &mut Vec<i32>,
    ) -> Result<IcrfKm, SpkKernelError> {
        if target_id == SSB_ID {
            return Ok(zero());
        }
        if visited.contains(&target_id) {
            return Err(SpkKernelError::CenterCycle { target_id });
        }
        visited.push(target_id);

        let mut out_of_range = None;
        for segment in self
            .segments
            .iter()
            .rev()
            .filter(|segment| segment.target_id == target_id)
        {
            if !segment.covers(jd) {
                out_of_range.get_or_insert(segment.out_of_range(jd));
                continue;
            }
            match segment.try_position(jd) {
                Ok(offset) => match self.try_position_from_ssb(segment.center_id, jd, visited) {
                    Ok(center) => {
                        visited.pop();
                        return Ok(center + offset);
                    }
                    Err(err) => {
                        out_of_range.get_or_insert(err);
                    }
                },
                Err(err @ EphemerisError::OutOfRange { .. }) => {
                    out_of_range.get_or_insert(SpkKernelError::Ephemeris(err));
                }
                Err(err) => {
                    visited.pop();
                    return Err(SpkKernelError::Ephemeris(err));
                }
            }
        }

        if target_id == EARTH_ID {
            if let Ok(earth) = self.try_earth_from_emb(jd, visited) {
                visited.pop();
                return Ok(earth);
            }
        }

        visited.pop();
        Err(out_of_range.unwrap_or(SpkKernelError::MissingSegment { target_id }))
    }

    fn try_earth_from_emb(
        &self,
        jd: JulianDate,
        visited: &mut Vec<i32>,
    ) -> Result<IcrfKm, SpkKernelError> {
        let emb = self.try_position_from_ssb(EMB_ID, jd, visited)?;
        let moon_off = self.try_edge_position(MOON_ID, EMB_ID, jd)?;
        Ok(emb - moon_off.scale(1.0 / EARTH_MOON_RATIO))
    }

    fn try_edge_position(
        &self,
        target_id: i32,
        center_id: i32,
        jd: JulianDate,
    ) -> Result<IcrfKm, SpkKernelError> {
        let mut out_of_range = None;
        for segment in self
            .segments
            .iter()
            .rev()
            .filter(|segment| segment.target_id == target_id && segment.center_id == center_id)
        {
            if !segment.covers(jd) {
                out_of_range.get_or_insert(segment.out_of_range(jd));
                continue;
            }
            match segment.try_position(jd) {
                Ok(position) => return Ok(position),
                Err(err @ EphemerisError::OutOfRange { .. }) => {
                    out_of_range.get_or_insert(SpkKernelError::Ephemeris(err));
                }
                Err(err) => return Err(SpkKernelError::Ephemeris(err)),
            }
        }
        Err(out_of_range.unwrap_or(SpkKernelError::MissingSegment { target_id }))
    }

    #[cfg(test)]
    fn from_indexed_segments(segments: Vec<IndexedSegmentData>) -> Self {
        Self {
            segments: segments
                .into_iter()
                .map(KernelSegment::from_indexed)
                .collect(),
        }
    }
}

fn zero() -> IcrfKm {
    Displacement::new(
        Kilometers::new(0.0),
        Kilometers::new(0.0),
        Kilometers::new(0.0),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::spk::SegmentData;

    fn indexed(target_id: i32, center_id: i32, x_km: f64) -> IndexedSegmentData {
        let intlen = 1000.0 * SECONDS_PER_DAY;
        IndexedSegmentData {
            target_id,
            center_id,
            frame_id: 1,
            start_et: 0.0,
            end_et: intlen,
            data: SegmentData {
                data_type: 2,
                init: 0.0,
                intlen,
                rsize: 8,
                ncoeff: 2,
                n_records: 1,
                records: vec![intlen / 2.0, intlen / 2.0, x_km, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
        }
    }

    fn jd_mid() -> JulianDate {
        JulianDate::new(2_451_545.0 + 500.0)
    }

    #[test]
    fn direct_target_observer_state_resolves() {
        let set =
            SpkKernelSet::from_indexed_segments(vec![indexed(5, 0, 500.0), indexed(399, 0, 100.0)]);
        let state = set.try_geometric_state(5, 399, jd_mid()).unwrap();
        assert!((state.x().value() - 400.0).abs() < 1.0e-9);
    }

    #[test]
    fn planet_center_chain_uses_satellite_offset() {
        let set = SpkKernelSet::from_indexed_segments(vec![
            indexed(5, 0, 500.0),
            indexed(599, 5, 7.0),
            indexed(399, 0, 100.0),
        ]);
        let state = set
            .try_major_planet_geocentric(MajorPlanet::Jupiter, PlanetPoint::Center, jd_mid())
            .unwrap();
        assert!((state.x().value() - 407.0).abs() < 1.0e-9);
    }

    #[test]
    fn missing_planet_center_is_explicit() {
        let set =
            SpkKernelSet::from_indexed_segments(vec![indexed(5, 0, 500.0), indexed(399, 0, 100.0)]);
        let err = set
            .try_major_planet_geocentric(MajorPlanet::Jupiter, PlanetPoint::Center, jd_mid())
            .unwrap_err();
        assert!(matches!(
            err,
            SpkKernelError::MissingSegment { target_id: 599 }
        ));
    }

    #[test]
    fn later_kernel_segment_wins() {
        let set = SpkKernelSet::from_indexed_segments(vec![
            indexed(5, 0, 500.0),
            indexed(5, 0, 900.0),
            indexed(399, 0, 100.0),
        ]);
        let state = set.try_geometric_state(5, 399, jd_mid()).unwrap();
        assert!((state.x().value() - 800.0).abs() < 1.0e-9);
    }

    #[test]
    fn out_of_coverage_epoch_is_explicit() {
        let set = SpkKernelSet::from_indexed_segments(vec![indexed(5, 0, 500.0)]);
        let err = set
            .try_geometric_state(5, 0, JulianDate::new(2_451_545.0 + 2000.0))
            .unwrap_err();
        assert!(matches!(
            err,
            SpkKernelError::Ephemeris(EphemerisError::OutOfRange { .. })
        ));
    }
}
