// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Sun-Earth Lagrange Chebyshev centers
//!
//! ## Scientific scope
//!
//! This module represents the instantaneous generalized Sun-Earth L1-L5
//! equilibrium points in a Sun-Earth rotating frame perturbed by the Moon. The
//! embedded archive is intentionally empty in this source tree until an official
//! DE440 generation run is committed; callers can still use the solver and fit
//! helpers to generate short fixtures or operational coefficient archives.
//!
//! ## Technical scope
//!
//! The public evaluator returns barycentric [`Position`] values in
//! [`AstronomicalUnit`] on the [`EclipticMeanJ2000`] frame from flat Chebyshev
//! records with layout `[mid_seconds, radius_seconds, x..., y..., z...]`, where
//! stored coefficients are kilometres. [`solver`] computes per-epoch positions
//! from an [`Ephemeris`] or [`DynEphemeris`], and [`fit`] builds records for a
//! caller-provided epoch span.
//!
//! ## References
//!
//! - Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*.
//! - Koon, W. S., Lo, M. W., Marsden, J. E., Ross, S. D. (2011). *Dynamical Systems,
//!   the Three-Body Problem and Space Mission Design*.
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440
//!   and DE441". *The Astronomical Journal* 161, 105.

pub mod fit;
pub mod solver;

use crate::coordinates::cartesian::Position;
use crate::coordinates::centers::Barycentric;
use crate::coordinates::frames::EclipticMeanJ2000;
use crate::embedded_data::lagrange as data;
use crate::ephemeris::EphemerisError;
use crate::qtty::{AstronomicalUnit, Days, Kilometer, Kilometers, Meters, Second};
use crate::time::JulianDate;

const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;
const J2000_JD: f64 = 2_451_545.0;

/// Sun-Earth Lagrange point selector.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SunEarthLagrangePoint {
    /// Collinear point between the Sun and Earth.
    L1,
    /// Collinear point beyond Earth, away from the Sun.
    L2,
    /// Collinear point beyond the Sun, opposite Earth.
    L3,
    /// Leading triangular point.
    L4,
    /// Trailing triangular point.
    L5,
}

impl SunEarthLagrangePoint {
    /// Returns the stable display label for this point.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::L1 => "L1",
            Self::L2 => "L2",
            Self::L3 => "L3",
            Self::L4 => "L4",
            Self::L5 => "L5",
        }
    }
}

/// Metadata describing the embedded Sun-Earth Lagrange Chebyshev archive.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LagrangeMetadata {
    /// Source ephemeris used by the generator.
    pub source: &'static str,
    /// First covered TT/TDB-like Julian Date.
    pub valid_from: Days,
    /// Last covered TT/TDB-like Julian Date.
    pub valid_to: Days,
    /// Name of the output reference frame.
    pub frame_name: &'static str,
    /// Name of the time scale used by the archive.
    pub time_scale_name: &'static str,
    /// Name of the length unit used by stored coefficients.
    pub length_unit_name: &'static str,
    /// Nominal Chebyshev block duration.
    pub block: Second,
    /// Validation grid step used by the generator.
    pub validation_step: Second,
    /// Maximum absolute validation error.
    pub max_abs_error: Meters,
    /// Generator version string.
    pub generator_version: &'static str,
    /// Generation timestamp.
    pub generated_at: &'static str,
    /// Archive checksum.
    pub checksum: &'static str,
}

/// Metadata for the embedded Sun-Earth Lagrange archive.
pub const SUN_EARTH_LAGRANGE_METADATA: LagrangeMetadata = LagrangeMetadata {
    source: "NONE - placeholder",
    valid_from: Days::new(J2000_JD),
    valid_to: Days::new(J2000_JD),
    frame_name: "EclipticMeanJ2000",
    time_scale_name: "TDB-compatible Julian Date",
    length_unit_name: "Kilometer",
    block: Second::new(0.0),
    validation_step: Second::new(0.0),
    max_abs_error: Meters::new(0.0),
    generator_version: "not generated",
    generated_at: "not generated",
    checksum: "empty",
};

/// Fallibly evaluates an embedded Sun-Earth Lagrange point archive record.
///
/// # Arguments
///
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
///
/// # Returns
///
/// Barycentric ecliptic J2000 position in astronomical units.
///
/// # Errors
///
/// Returns [`EphemerisError::OutOfRange`] when the embedded archive is empty or
/// the requested epoch lies outside its records.
pub fn try_sun_earth_lagrange_barycentric(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    evaluate_embedded(point, jd)
}

/// Evaluates an embedded Sun-Earth Lagrange point archive record.
///
/// # Arguments
///
/// - `point`: L1-L5 selector.
/// - `jd`: Epoch as a typed Julian Date.
///
/// # Returns
///
/// Barycentric ecliptic J2000 position in astronomical units.
///
/// # Panics
///
/// Panics when the embedded archive is empty or the requested epoch lies outside
/// coverage, matching the infallible JPL backend accessors.
pub fn sun_earth_lagrange_barycentric(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_sun_earth_lagrange_barycentric(point, jd)
        .expect("Sun-Earth Lagrange archive epoch outside coverage")
}

pub(crate) fn evaluate_records(
    records: &[f64],
    ncoeff: usize,
    jd: JulianDate,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let record_len = 2 + 3 * ncoeff;
    let jd_value = jd.raw().value();
    if records.is_empty()
        || ncoeff == 0
        || record_len == 2
        || !records.len().is_multiple_of(record_len)
    {
        return Err(EphemerisError::OutOfRange {
            jd: jd_value,
            start_jd: SUN_EARTH_LAGRANGE_METADATA.valid_from.value(),
            end_jd: SUN_EARTH_LAGRANGE_METADATA.valid_to.value(),
        });
    }

    let seconds = (jd_value - J2000_JD) * SECONDS_PER_DAY;
    let mut first = f64::INFINITY;
    let mut last = f64::NEG_INFINITY;
    for record in records.chunks_exact(record_len) {
        let start = record[0] - record[1];
        let end = record[0] + record[1];
        first = first.min(start);
        last = last.max(end);
        if seconds >= start && seconds <= end {
            let tau = (seconds - record[0]) / record[1];
            let cx = &record[2..2 + ncoeff];
            let cy = &record[2 + ncoeff..2 + 2 * ncoeff];
            let cz = &record[2 + 2 * ncoeff..2 + 3 * ncoeff];
            let km = Position::<Barycentric, EclipticMeanJ2000, Kilometer>::new(
                Kilometers::new(cheby::evaluate(cx, tau)),
                Kilometers::new(cheby::evaluate(cy, tau)),
                Kilometers::new(cheby::evaluate(cz, tau)),
            );
            return Ok(km.to_unit::<AstronomicalUnit>());
        }
    }

    Err(EphemerisError::OutOfRange {
        jd: jd_value,
        start_jd: J2000_JD + first / SECONDS_PER_DAY,
        end_jd: J2000_JD + last / SECONDS_PER_DAY,
    })
}

fn evaluate_embedded(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    evaluate_records(records_for(point), data::NCOEFF, jd)
}

fn records_for(point: SunEarthLagrangePoint) -> &'static [f64] {
    match point {
        SunEarthLagrangePoint::L1 => data::RECORDS_L1,
        SunEarthLagrangePoint::L2 => data::RECORDS_L2,
        SunEarthLagrangePoint::L3 => data::RECORDS_L3,
        SunEarthLagrangePoint::L4 => data::RECORDS_L4,
        SunEarthLagrangePoint::L5 => data::RECORDS_L5,
    }
}

pub(crate) fn fallback_or_solve<Eph: crate::ephemeris::Ephemeris>(
    point: SunEarthLagrangePoint,
    jd: JulianDate,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_sun_earth_lagrange_barycentric(point, jd).unwrap_or_else(|_| {
        solver::solve_sun_earth_lagrange::<Eph>(point, jd)
            .expect("Sun-Earth Lagrange solver failed")
            .position
            .to_unit::<AstronomicalUnit>()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::J2000;

    #[test]
    fn empty_archive_is_out_of_range() {
        let err = try_sun_earth_lagrange_barycentric(SunEarthLagrangePoint::L1, J2000)
            .expect_err("placeholder archive must not evaluate");
        assert!(matches!(err, EphemerisError::OutOfRange { .. }));
    }

    #[test]
    fn point_labels_are_stable() {
        assert_eq!(SunEarthLagrangePoint::L5.label(), "L5");
    }

    #[test]
    fn astronomical_units_constructor_is_available() {
        let zero = crate::qtty::AstronomicalUnits::new(0.0);
        assert_eq!(zero.value(), 0.0);
    }
}
