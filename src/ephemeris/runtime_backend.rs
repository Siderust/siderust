// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime-loaded ephemeris backend.
//!
//! [`RuntimeEphemeris`] loads a JPL DE4xx BSP file at runtime and provides
//! the same body-chain computations as the compile-time backends.
//! It implements [`DynEphemeris`](super::DynEphemeris) (instance-based,
//! object-safe).

use super::{AuPerDay, DynEphemeris, EphemerisError};
use crate::archive::ArchiveError;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::ephemeris::jpl::bodies;
use crate::ephemeris::jpl::eval::DynSegmentStack;
use crate::formats::spice::{self, spk};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::JulianDate;
use std::path::Path;
use std::sync::Arc;

/// Shared inner data for a runtime-loaded ephemeris.
struct RuntimeEphemerisInner {
    sun: DynSegmentStack,
    emb: DynSegmentStack,
    moon: DynSegmentStack,
    earth: Option<DynSegmentStack>,
}

/// Runtime-loaded JPL DE4xx ephemeris backend.
///
/// This struct loads a BSP file at runtime and evaluates Sun, Earth, and Moon
/// positions using the same Chebyshev polynomial evaluation as the compile-time
/// backends. Unlike compile-time VSOP87, it:
///
/// - Does **not** require a Cargo feature flag
/// - Stores coefficient data on the **heap** (via `Vec<f64>`)
/// - Implements [`DynEphemeris`] (instance methods with `&self`)
/// - Is cloneable and shareable via internal `Arc`
/// - Retains **all** matching SPK segments per body chain and selects by epoch
///   (later kernels win when coverage overlaps), matching [`crate::formats::spice::SpkKernelSet`] semantics
///
/// # Example
///
/// ```rust,ignore
/// use siderust::ephemeris::{RuntimeEphemeris, DynEphemeris};
/// use siderust::time::JulianDate;
///
/// let eph = RuntimeEphemeris::from_bsp("path/to/de440.bsp")?;
/// let sun_pos = eph.sun_barycentric(crate::J2000);
/// ```
#[derive(Clone)]
pub struct RuntimeEphemeris {
    inner: Arc<RuntimeEphemerisInner>,
}

impl RuntimeEphemeris {
    /// Load a runtime ephemeris from a BSP file on disk.
    ///
    /// The file is read entirely into memory, parsed as a DAF/SPK container,
    /// and every supported Sun, EMB, Moon, and optional Earth segment is indexed.
    pub fn from_bsp(path: impl AsRef<Path>) -> Result<Self, ArchiveError> {
        let file_data = std::fs::read(path.as_ref())?;
        Self::from_bytes(&file_data)
    }

    /// Load a runtime ephemeris from raw BSP bytes already in memory.
    pub fn from_bytes(data: &[u8]) -> Result<Self, ArchiveError> {
        let indexed = spk::parse_indexed_segments(data).map_err(spice_error_to_archive)?;
        Ok(Self::from_indexed_segments(indexed))
    }

    /// Construct from every parsed J2000 Type 2/3 segment in a BSP file.
    pub fn from_indexed_segments(indexed: Vec<spk::IndexedSegmentData>) -> Self {
        let inner = RuntimeEphemerisInner {
            sun: DynSegmentStack::for_indexed_pair(&indexed, spk::SUN_TARGET, spk::SUN_CENTER),
            emb: DynSegmentStack::for_indexed_pair(&indexed, spk::EMB_TARGET, spk::EMB_CENTER),
            moon: DynSegmentStack::for_indexed_pair(&indexed, spk::MOON_TARGET, spk::MOON_CENTER),
            earth: {
                let stack = DynSegmentStack::for_indexed_pair(
                    &indexed,
                    spk::EARTH_TARGET,
                    spk::EARTH_CENTER,
                );
                if stack.segment_count() > 0 {
                    Some(stack)
                } else {
                    None
                }
            },
        };
        Self {
            inner: Arc::new(inner),
        }
    }

    /// Construct from legacy single-segment [`spk::BspSegments`] (tests and tooling).
    pub fn from_segments(segments: spk::BspSegments) -> Self {
        let inner = RuntimeEphemerisInner {
            sun: DynSegmentStack::from_spk_segment(&segments.sun),
            emb: DynSegmentStack::from_spk_segment(&segments.emb),
            moon: DynSegmentStack::from_spk_segment(&segments.moon),
            earth: segments
                .earth
                .as_ref()
                .map(DynSegmentStack::from_spk_segment),
        };
        Self {
            inner: Arc::new(inner),
        }
    }

    /// Load from a BSP file resolved by [`siderust_archive::jpl::DatasetManager`].
    ///
    /// ```rust,ignore
    /// use siderust_archive::jpl::{DatasetManager, refs::JplDatasetId};
    /// use siderust::ephemeris::RuntimeEphemeris;
    ///
    /// let dm = DatasetManager::new()?;
    /// let eph = RuntimeEphemeris::from_dataset_manager(&dm, JplDatasetId::De441)?;
    /// ```
    #[cfg(feature = "runtime-data")]
    pub fn from_dataset_manager(
        dm: &crate::archive::jpl::DatasetManager,
        id: crate::archive::jpl::JplDatasetId,
    ) -> Result<Self, ArchiveError> {
        let path = dm.ensure(id)?;
        Self::from_bsp(path)
    }
}

fn spice_error_to_archive(err: spice::SpiceError) -> ArchiveError {
    match err {
        spice::SpiceError::Io(e) => ArchiveError::Io(e),
        other => ArchiveError::Integrity(format!("SPICE parse error: {other}")),
    }
}

impl std::fmt::Debug for RuntimeEphemeris {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RuntimeEphemeris")
            .field("sun_segments", &self.inner.sun.segment_count())
            .field("emb_segments", &self.inner.emb.segment_count())
            .field("moon_segments", &self.inner.moon.segment_count())
            .field(
                "earth_segments",
                &self
                    .inner
                    .earth
                    .as_ref()
                    .map(DynSegmentStack::segment_count),
            )
            .finish()
    }
}

impl DynEphemeris for RuntimeEphemeris {
    #[inline]
    fn try_sun_barycentric(
        &self,
        jd: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        bodies::try_dyn_sun_barycentric(jd, &self.inner.sun)
    }

    #[inline]
    fn sun_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::dyn_sun_barycentric(jd, &self.inner.sun)
    }

    #[inline]
    fn try_earth_barycentric(
        &self,
        jd: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_barycentric_direct(jd, &self.inner.emb, earth);
        }
        bodies::try_dyn_earth_barycentric(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn earth_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_barycentric_direct(jd, &self.inner.emb, earth)
                .expect("runtime JPL Earth barycentric position unavailable");
        }
        bodies::dyn_earth_barycentric(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn try_earth_heliocentric(
        &self,
        jd: JulianDate,
    ) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_heliocentric_direct(
                jd,
                &self.inner.sun,
                &self.inner.emb,
                earth,
            );
        }
        bodies::try_dyn_earth_heliocentric(jd, &self.inner.sun, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn earth_heliocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_heliocentric_direct(
                jd,
                &self.inner.sun,
                &self.inner.emb,
                earth,
            )
            .expect("runtime JPL Earth heliocentric position unavailable");
        }
        bodies::dyn_earth_heliocentric(jd, &self.inner.sun, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn try_earth_barycentric_velocity(
        &self,
        jd: JulianDate,
    ) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_barycentric_velocity_direct(jd, &self.inner.emb, earth);
        }
        bodies::try_dyn_earth_barycentric_velocity(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn earth_barycentric_velocity(&self, jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_earth_barycentric_velocity_direct(jd, &self.inner.emb, earth)
                .expect("runtime JPL Earth barycentric velocity unavailable");
        }
        bodies::dyn_earth_barycentric_velocity(jd, &self.inner.emb, &self.inner.moon)
    }

    #[inline]
    fn try_moon_geocentric(
        &self,
        jd: JulianDate,
    ) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_moon_geocentric_direct(jd, &self.inner.moon, earth);
        }
        bodies::try_dyn_moon_geocentric(jd, &self.inner.moon)
    }

    #[inline]
    fn moon_geocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        if let Some(earth) = &self.inner.earth {
            return bodies::try_dyn_moon_geocentric_direct(jd, &self.inner.moon, earth)
                .expect("runtime JPL Moon geocentric position unavailable");
        }
        bodies::dyn_moon_geocentric(jd, &self.inner.moon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::spice::spk::{BspSegments, IndexedSegmentData, SegmentData};

    const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;
    const JD_J2000: f64 = tempoch::J2000_JD_TT_DAY.value();

    /// Create a minimal SegmentData with constant position (x_km, y_km, z_km).
    fn make_segment(x_km: f64, y_km: f64, z_km: f64) -> SegmentData {
        let ncoeff = 2usize;
        let rsize = 2 + 3 * ncoeff;
        let intlen = 1000.0 * SECONDS_PER_DAY;
        let mid = intlen / 2.0;
        let radius = intlen / 2.0;
        let records = vec![mid, radius, x_km, 0.0, y_km, 0.0, z_km, 0.0];
        SegmentData {
            data_type: 2,
            init: 0.0,
            intlen,
            rsize,
            ncoeff,
            n_records: 1,
            records,
        }
    }

    fn indexed(
        target_id: i32,
        center_id: i32,
        x_km: f64,
        start_et: f64,
        end_et: f64,
    ) -> IndexedSegmentData {
        let intlen = end_et - start_et;
        IndexedSegmentData {
            target_id,
            center_id,
            frame_id: 1,
            start_et,
            end_et,
            data: {
                let mut data = make_segment(x_km, 0.0, 0.0);
                data.init = start_et;
                data.intlen = intlen;
                data.records[0] = start_et + intlen / 2.0;
                data.records[1] = intlen / 2.0;
                data
            },
        }
    }

    fn make_bsp_segments() -> BspSegments {
        BspSegments {
            sun: make_segment(1.0e8, 2.0e7, 1.0e6),
            emb: make_segment(1.5e8, 0.0, 0.0),
            moon: make_segment(3.84e5, 5.0e3, 1.0e3),
            earth: None,
        }
    }

    fn make_bsp_segments_with_earth() -> BspSegments {
        BspSegments {
            sun: make_segment(1.0e8, 0.0, 0.0),
            emb: make_segment(1.5e8, 0.0, 0.0),
            moon: make_segment(3.84e5, 0.0, 0.0),
            earth: Some(make_segment(-5.0e3, 0.0, 0.0)),
        }
    }

    fn jd_mid() -> JulianDate {
        crate::time::JulianDate::new(JD_J2000 + 500.0)
    }

    #[test]
    fn from_segments_roundtrip_n_records() {
        let segs = make_bsp_segments();
        let eph = RuntimeEphemeris::from_segments(segs);
        let dbg = format!("{eph:?}");
        assert!(dbg.contains("RuntimeEphemeris"));
        assert!(dbg.contains("sun_segments"));
    }

    #[test]
    fn clone_gives_same_results() {
        let segs = make_bsp_segments();
        let eph = RuntimeEphemeris::from_segments(segs);
        let eph2 = eph.clone();
        let jd = jd_mid();
        let pos1 = eph.sun_barycentric(jd);
        let pos2 = eph2.sun_barycentric(jd);
        assert!((pos1.x().value() - pos2.x().value()).abs() < 1e-15);
    }

    #[test]
    fn sun_barycentric_is_finite() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments());
        let pos = eph.sun_barycentric(jd_mid());
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_barycentric_is_finite() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments());
        let pos = eph.earth_barycentric(jd_mid());
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_heliocentric_is_finite() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments());
        let pos = eph.earth_heliocentric(jd_mid());
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_barycentric_velocity_is_finite() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments());
        let vel = eph.earth_barycentric_velocity(jd_mid());
        assert!(vel.x().is_finite());
        assert!(vel.y().is_finite());
        assert!(vel.z().is_finite());
    }

    #[test]
    fn moon_geocentric_is_finite() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments());
        let pos = eph.moon_geocentric(jd_mid());
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_segment_takes_precedence_for_earth_barycentric() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments_with_earth());
        let pos = eph.earth_barycentric(jd_mid());
        let mag_au =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected_au = (1.5e8 - 5.0e3) / 149_597_870.700;

        assert!(
            (mag_au - expected_au).abs() < 1e-12,
            "mag_au={mag_au}, expected={expected_au}"
        );
    }

    #[test]
    fn earth_segment_takes_precedence_for_moon_geocentric() {
        let eph = RuntimeEphemeris::from_segments(make_bsp_segments_with_earth());
        let pos = eph.moon_geocentric(jd_mid());
        let mag_km =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected_km = 384_000.0 + 5_000.0;

        assert!(
            (mag_km - expected_km).abs() < 1e-8,
            "mag_km={mag_km}, expected={expected_km}"
        );
    }

    #[test]
    fn later_sun_segment_wins_at_overlap() {
        let intlen = 1000.0 * SECONDS_PER_DAY;
        let indexed = vec![
            indexed(spk::SUN_TARGET, spk::SUN_CENTER, 100.0, 0.0, intlen),
            indexed(spk::SUN_TARGET, spk::SUN_CENTER, 900.0, 0.0, intlen),
            indexed(spk::EMB_TARGET, spk::EMB_CENTER, 1.5e8, 0.0, intlen),
            indexed(spk::MOON_TARGET, spk::MOON_CENTER, 3.84e5, 0.0, intlen),
        ];
        let eph = RuntimeEphemeris::from_indexed_segments(indexed);
        let pos = eph.sun_barycentric(jd_mid());
        let mag_km = (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2))
            .sqrt()
            * 149_597_870.700;
        assert!(
            (mag_km - 900.0).abs() < 1.0,
            "expected later Sun segment (900 km), got {mag_km} km"
        );
    }

    #[test]
    fn direct_earth_matches_emb_moon_derivation() {
        let intlen = 1000.0 * SECONDS_PER_DAY;
        let earth_off_km = -4.0e3;
        let moon_off_km = -earth_off_km * crate::archive::jpl::constants::EARTH_MOON_RATIO;
        let with_earth = RuntimeEphemeris::from_indexed_segments(vec![
            indexed(spk::SUN_TARGET, spk::SUN_CENTER, 1.0e8, 0.0, intlen),
            indexed(spk::EMB_TARGET, spk::EMB_CENTER, 1.5e8, 0.0, intlen),
            indexed(spk::MOON_TARGET, spk::MOON_CENTER, moon_off_km, 0.0, intlen),
            indexed(
                spk::EARTH_TARGET,
                spk::EARTH_CENTER,
                earth_off_km,
                0.0,
                intlen,
            ),
        ]);
        let without_earth = RuntimeEphemeris::from_indexed_segments(vec![
            indexed(spk::SUN_TARGET, spk::SUN_CENTER, 1.0e8, 0.0, intlen),
            indexed(spk::EMB_TARGET, spk::EMB_CENTER, 1.5e8, 0.0, intlen),
            indexed(spk::MOON_TARGET, spk::MOON_CENTER, moon_off_km, 0.0, intlen),
        ]);
        let jd = jd_mid();
        let direct = with_earth.earth_barycentric(jd);
        let derived = without_earth.earth_barycentric(jd);
        assert!((direct.x().value() - derived.x().value()).abs() < 1e-9);
    }

    #[test]
    fn from_bytes_on_invalid_data_returns_error() {
        let data = b"not a bsp file";
        let result = RuntimeEphemeris::from_bytes(data);
        assert!(result.is_err());
    }

    #[test]
    fn from_bsp_on_nonexistent_file_returns_error() {
        let result = RuntimeEphemeris::from_bsp("/nonexistent/path/de999.bsp");
        assert!(result.is_err());
    }
}
