// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime-loaded ephemeris backend.
//!
//! [`RuntimeEphemeris`] loads a JPL DE4xx BSP file at runtime and provides
//! the same body-chain computations as the compile-time backends.
//! It implements [`DynEphemeris`](super::DynEphemeris) (instance-based,
//! object-safe).

use super::{AuPerDay, DynEphemeris};
use crate::calculus::jpl::bodies;
use crate::calculus::jpl::eval::DynSegmentDescriptor;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::data::spk;
use crate::data::DataError;
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, Kilometer};
use std::path::Path;
use std::sync::Arc;

/// Shared inner data for a runtime-loaded ephemeris.
struct RuntimeEphemerisInner {
    sun: DynSegmentDescriptor,
    emb: DynSegmentDescriptor,
    moon: DynSegmentDescriptor,
}

/// Runtime-loaded JPL DE4xx ephemeris backend.
///
/// This struct loads a BSP file at runtime and evaluates Sun, Earth, and Moon
/// positions using the same Chebyshev polynomial evaluation as the compile-time
/// backends. Unlike `De440Ephemeris`/`De441Ephemeris`, it:
///
/// - Does **not** require a Cargo feature flag
/// - Stores coefficient data on the **heap** (via `Vec<f64>`)
/// - Implements [`DynEphemeris`] (instance methods with `&self`)
/// - Is cloneable and shareable via internal `Arc`
///
/// # Example
///
/// ```rust,ignore
/// use siderust::calculus::ephemeris::{RuntimeEphemeris, DynEphemeris};
/// use siderust::time::JulianDate;
///
/// let eph = RuntimeEphemeris::from_bsp("path/to/de440.bsp")?;
/// let sun_pos = eph.sun_barycentric(JulianDate::J2000);
/// ```
#[derive(Clone)]
pub struct RuntimeEphemeris {
    inner: Arc<RuntimeEphemerisInner>,
}

impl RuntimeEphemeris {
    /// Load a runtime ephemeris from a BSP file on disk.
    ///
    /// The file is read entirely into memory, parsed as a DAF/SPK container,
    /// and the Sun, EMB, and Moon Chebyshev segments are extracted.
    pub fn from_bsp(path: impl AsRef<Path>) -> Result<Self, DataError> {
        let file_data = std::fs::read(path.as_ref())?;
        Self::from_bytes(&file_data)
    }

    /// Load a runtime ephemeris from raw BSP bytes already in memory.
    pub fn from_bytes(data: &[u8]) -> Result<Self, DataError> {
        let segments = spk::parse_bsp(data)?;
        Ok(Self::from_segments(segments))
    }

    /// Construct from pre-parsed BSP segments.
    pub fn from_segments(segments: spk::BspSegments) -> Self {
        let inner = RuntimeEphemerisInner {
            sun: DynSegmentDescriptor::from_spk(&segments.sun),
            emb: DynSegmentDescriptor::from_spk(&segments.emb),
            moon: DynSegmentDescriptor::from_spk(&segments.moon),
        };
        Self {
            inner: Arc::new(inner),
        }
    }

    /// Load from a BSP file resolved by the [`DataManager`](crate::data::DataManager).
    ///
    /// This is a convenience method that combines data management with loading:
    /// ```rust,ignore
    /// use siderust::data::{DataManager, DatasetId};
    /// use siderust::calculus::ephemeris::RuntimeEphemeris;
    ///
    /// let dm = DataManager::new()?;
    /// let path = dm.ensure(DatasetId::De441)?;
    /// let eph = RuntimeEphemeris::from_bsp(path)?;
    /// ```
    #[cfg(feature = "runtime-data")]
    pub fn from_data_manager(
        dm: &crate::data::DataManager,
        id: crate::data::DatasetId,
    ) -> Result<Self, DataError> {
        let path = dm.ensure(id)?;
        Self::from_bsp(path)
    }
}

impl std::fmt::Debug for RuntimeEphemeris {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RuntimeEphemeris")
            .field("sun_records", &self.inner.sun.n_records)
            .field("emb_records", &self.inner.emb.n_records)
            .field("moon_records", &self.inner.moon.n_records)
            .finish()
    }
}

impl DynEphemeris for RuntimeEphemeris {
    #[inline]
    fn sun_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::dyn_sun_barycentric(jd, &self.inner.sun)
    }

    #[inline]
    fn earth_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::dyn_earth_barycentric(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn earth_heliocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::dyn_earth_heliocentric(jd, &self.inner.sun, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn earth_barycentric_velocity(&self, jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        bodies::dyn_earth_barycentric_velocity(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn moon_geocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        bodies::dyn_moon_geocentric(jd, &self.inner.moon)
    }
}
