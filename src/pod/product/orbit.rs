//! # Orbit product assembly
//!
//! ## Scientific scope
//!
//! Precise orbit determination ultimately produces a time history of
//! spacecraft states that must be exchanged in community-standard orbit
//! formats. This module bridges that scientific result into SP3 and OEM
//! products without changing the underlying trajectory meaning.
//!
//! The validity of the output therefore depends entirely on the supplied
//! state series and metadata. The module does not smooth, interpolate, or
//! otherwise alter the orbit solution.
//!
//! ## Technical scope
//!
//! The public entry points are `write_sp3_from_states` and
//! `write_oem_from_states`. They accept a sequence of `OrbitState` values
//! and the minimal metadata needed by the target file format, then delegate
//! low-level serialization to `siderust::formats`.
//!
//! A convenience wrapper `write_oem_from_spacecraft_states` accepts
//! `SpacecraftState` values directly by extracting the inner `orbit` field.
//!
//! Residual statistics, manifests, and QC summaries are intentionally
//! handled elsewhere.
//!
//! ## References
//!
//! - Consultative Committee for Space Data Systems. (2010). Orbit Data
//!   Messages, CCSDS 502.0-B-2 / 502.0-B-3.
//! - International GNSS Service. (2020). SP3-c / SP3-d Orbit Format
//!   Specification.
use crate::astro::dynamics::state::SpacecraftState;
use crate::astro::dynamics::OrbitState;
use crate::coordinates::frames::GCRS;
use crate::formats::ccsds::oem::{write_oem, OemMetadata, OemState};
use crate::formats::igs::sp3::{write_sp3, EarthCenter, Sp3Epoch, Sp3Position, Sp3Record};
use crate::formats::FormatError;
use crate::time::JulianDate;
use affn::cartesian;
use qtty::time::Microseconds;
use qtty::unit::Kilometer;
use qtty::Day;
use std::io::Write;
use std::path::Path;
use tempoch::{Time, UTC};

use super::error::PodProductsError;

/// Build an SP3 record from an in-memory state series and write it to `w`.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::orbit::write_sp3_from_states;
/// use siderust::astro::dynamics::{OrbitState, Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// let states: Vec<OrbitState> = (0..3).map(|i| {
///     OrbitState::new(
///         JulianDate::new(2_451_545.0 + i as f64 * 30.0 / 86_400.0).to_j2000s(),
///         Position::<GCRS>::new(7000.0, 0.0, 0.0),
///         Velocity::<GCRS>::new(0.0, 7.5, 0.0),
///     )
/// }).collect();
/// let mut buf = Vec::new();
/// write_sp3_from_states(&mut buf, "L01", &states).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert_eq!(text.lines().filter(|l| l.starts_with("PL01")).count(), 3);
/// ```
pub fn write_sp3_from_states<W: Write>(
    w: &mut W,
    sat_id: &str,
    states: &[OrbitState],
) -> Result<(), crate::formats::igs::sp3::Sp3Error> {
    let header = vec![
        format!(
            "#dP2024  1  1  0  0  0.00000000      {:5} ORBIT IGS20 HLM  POD",
            states.len()
        ),
        "##  2295      0.00000000   900.00000000 60310 0.0000000000000".to_string(),
        format!("+    1   {sat_id:<3}                                                     "),
        "++         5                                                                  "
            .to_string(),
        "%c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc".to_string(),
        "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc".to_string(),
        "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000".to_string(),
        "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000".to_string(),
        "%i    0    0    0    0      0      0      0      0         0".to_string(),
        "%i    0    0    0    0      0      0      0      0         0".to_string(),
        "/* siderust-pod orbit product".to_string(),
    ];
    let epochs = states
        .iter()
        .map(|s| {
            // SP3 epochs are UTC; orbit states carry TT internally.
            let epoch_utc: Time<UTC> =
                JulianDate::try_new(Day::new(s.epoch.to::<tempoch::JD>().value()))
                    .expect("OrbitState epoch must be finite")
                    .to_scale::<UTC>()
                    .into();
            Sp3Epoch {
                time: epoch_utc,
                positions: vec![Sp3Position {
                    sat_id: sat_id.to_string(),
                    // Sp3Position.position is cartesian::Position<EarthCenter, GCRS, Kilometer>
                    // where EarthCenter is crate::formats::igs::sp3::EarthCenter (Params = ()).
                    position: cartesian::Position::<EarthCenter, GCRS, Kilometer>::new(
                        s.position.x().value(),
                        s.position.y().value(),
                        s.position.z().value(),
                    ),
                    clock: Microseconds::new(999_999.999_999),
                }],
            }
        })
        .collect();
    let rec = Sp3Record { header, epochs };
    write_sp3(w, &rec)
}

/// Write an OEM product from an `OrbitState` series.
///
/// Converts each [`OrbitState`] into an [`OemState`] (km/km·s⁻¹, TT Julian
/// Date) and delegates to `siderust::formats::ccsds::oem::write_oem`.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::orbit::write_oem_from_states;
/// use siderust::astro::dynamics::{OrbitState, Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// let states: Vec<OrbitState> = (0..3).map(|i| {
///     OrbitState::new(
///         JulianDate::new(2_451_545.0 + i as f64 * 30.0 / 86_400.0).to_j2000s(),
///         Position::<GCRS>::new(7000.0, 0.0, 0.0),
///         Velocity::<GCRS>::new(0.0, 7.5, 0.0),
///     )
/// }).collect();
/// let mut buf = Vec::new();
/// write_oem_from_states(&mut buf, "1900-001A", "TEST", &states).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert_eq!(text.lines().filter(|l| l.starts_with("2000-")).count(), 3);
/// ```
pub fn write_oem_from_states<W: Write>(
    w: &mut W,
    object_id: &str,
    object_name: &str,
    states: &[OrbitState],
) -> Result<(), FormatError> {
    let meta = OemMetadata {
        object_id: object_id.into(),
        object_name: object_name.into(),
        ref_frame: "EME2000".into(),
        time_system: "TT".into(),
        center_name: "EARTH".into(),
    };
    let oem_states: Vec<OemState> = states
        .iter()
        .map(|s| {
            OemState::new(
                s.epoch.to::<tempoch::JD>().value(),
                [
                    s.position.x().value(),
                    s.position.y().value(),
                    s.position.z().value(),
                ],
                [
                    s.velocity.x().value(),
                    s.velocity.y().value(),
                    s.velocity.z().value(),
                ],
            )
        })
        .collect();
    write_oem(w, &meta, &oem_states)
}

/// Write an OEM product from a `SpacecraftState` series.
///
/// Convenience wrapper that extracts the inner `orbit` field from each
/// [`SpacecraftState`] and delegates to [`write_oem_from_states`].
///
/// # Examples
///
/// ```
/// use siderust::pod::product::orbit::write_oem_from_spacecraft_states;
/// use siderust::astro::dynamics::state::{SpacecraftProperties, SpacecraftState};
/// use siderust::astro::dynamics::{OrbitState, Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// let orbit = OrbitState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
/// );
/// let sc = SpacecraftState { orbit, properties: SpacecraftProperties::demo_leo() };
/// let mut buf = Vec::new();
/// write_oem_from_spacecraft_states(&mut buf, "1900-001A", "TEST", &[sc]).unwrap();
/// assert!(String::from_utf8(buf).unwrap().contains("OBJECT_ID"));
/// ```
pub fn write_oem_from_spacecraft_states<W: Write>(
    w: &mut W,
    object_id: &str,
    object_name: &str,
    states: &[SpacecraftState],
) -> Result<(), FormatError> {
    let orbits: Vec<OrbitState> = states.iter().map(|s| s.orbit).collect();
    write_oem_from_states(w, object_id, object_name, &orbits)
}

// ─── Writer wrappers ─────────────────────────────────────────────────────────

/// Wrapper around the SP3 writer that writes to a file path.
///
/// Carries the satellite identifier so the same writer instance can be
/// reused for multiple output paths. To write to an existing
/// [`std::io::Write`] sink use [`write_sp3_from_states`] directly.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::orbit::Sp3ProductWriter;
///
/// let w = Sp3ProductWriter::new("G01");
/// drop(w);
/// ```
pub struct Sp3ProductWriter {
    sat_id: String,
}

impl Sp3ProductWriter {
    /// Create a new writer for the given satellite identifier.
    ///
    /// `sat_id` is the 3-character SP3 satellite identifier (e.g. `"L01"`).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::orbit::Sp3ProductWriter;
    ///
    /// let w = Sp3ProductWriter::new("G01");
    /// drop(w);
    /// ```
    pub fn new(sat_id: impl Into<String>) -> Self {
        Self {
            sat_id: sat_id.into(),
        }
    }

    /// Write an SP3 product to `path` from the given `SpacecraftState` slice.
    ///
    /// Creates or truncates the file at `path`, then serialises all states.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Io`] if the file cannot be created or
    /// written, or [`PodProductsError::Sp3`] if SP3 serialisation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use siderust::pod::product::orbit::Sp3ProductWriter;
    /// use std::path::Path;
    ///
    /// // Given `states: Vec<SpacecraftState>`:
    /// // Sp3ProductWriter::new("L01")
    /// //     .write_from_states(Path::new("orbit.sp3"), &states)
    /// //     .unwrap();
    /// ```
    pub fn write_from_states(
        &self,
        path: &Path,
        states: &[SpacecraftState],
    ) -> Result<(), PodProductsError> {
        let orbits: Vec<OrbitState> = states.iter().map(|s| s.orbit).collect();
        let mut f = std::fs::File::create(path)?;
        write_sp3_from_states(&mut f, &self.sat_id, &orbits)?;
        Ok(())
    }
}

/// Wrapper around the OEM writer that writes to a file path.
///
/// Carries the CCSDS object metadata so the same writer instance can be
/// reused for multiple output paths. To write to an existing
/// [`std::io::Write`] sink use [`write_oem_from_spacecraft_states`] directly.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::orbit::OemProductWriter;
///
/// let w = OemProductWriter::new("2024-001A", "MYSAT");
/// drop(w);
/// ```
pub struct OemProductWriter {
    object_id: String,
    object_name: String,
}

impl OemProductWriter {
    /// Create a new writer with the given CCSDS object identifiers.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::orbit::OemProductWriter;
    ///
    /// let w = OemProductWriter::new("2024-001A", "MYSAT");
    /// drop(w);
    /// ```
    pub fn new(object_id: impl Into<String>, object_name: impl Into<String>) -> Self {
        Self {
            object_id: object_id.into(),
            object_name: object_name.into(),
        }
    }

    /// Write an OEM product to `path` from the given `SpacecraftState` slice.
    ///
    /// Creates or truncates the file at `path`, then serialises all states.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Io`] if the file cannot be created or
    /// written, or [`PodProductsError::PodIo`] if OEM serialisation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use siderust::pod::product::orbit::OemProductWriter;
    /// use std::path::Path;
    ///
    /// // Given `states: Vec<SpacecraftState>`:
    /// // OemProductWriter::new("2024-001A", "MYSAT")
    /// //     .write_from_states(Path::new("orbit.oem"), &states)
    /// //     .unwrap();
    /// ```
    pub fn write_from_states(
        &self,
        path: &Path,
        states: &[SpacecraftState],
    ) -> Result<(), PodProductsError> {
        let mut f = std::fs::File::create(path)?;
        write_oem_from_spacecraft_states(&mut f, &self.object_id, &self.object_name, states)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::state::{Position, SpacecraftProperties, Velocity};
    use crate::time::JulianDate;

    fn fake_orbit_states() -> Vec<OrbitState> {
        (0..5)
            .map(|i| {
                OrbitState::new(
                    JulianDate::new(2_451_545.0 + i as f64 * 30.0 / 86_400.0).to_j2000s(),
                    Position::<GCRS>::new(7000.0 + i as f64, 0.0, 0.0),
                    Velocity::<GCRS>::new(0.0, 7.5, 0.0),
                )
            })
            .collect()
    }

    fn fake_spacecraft_states() -> Vec<SpacecraftState> {
        fake_orbit_states()
            .into_iter()
            .map(|orbit| SpacecraftState {
                orbit,
                properties: SpacecraftProperties::demo_leo(),
            })
            .collect()
    }

    #[test]
    fn sp3_writer_emits_n_epochs() {
        let s = fake_orbit_states();
        let mut buf = Vec::new();
        write_sp3_from_states(&mut buf, "L01", &s).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with("PL01")).count(), 5);
    }

    #[test]
    fn oem_writer_emits_data_lines() {
        let s = fake_orbit_states();
        let mut buf = Vec::new();
        write_oem_from_states(&mut buf, "1900-001A", "TEST", &s).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with("2000-")).count(), 5);
    }

    #[test]
    fn oem_from_spacecraft_states_writes_object_id() {
        let scs = fake_spacecraft_states();
        let mut buf = Vec::new();
        write_oem_from_spacecraft_states(&mut buf, "2024-001A", "MYSAT", &scs).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("OBJECT_ID            = 2024-001A"));
        assert!(text.contains("OBJECT_NAME          = MYSAT"));
    }

    #[test]
    fn sp3_product_writer_roundtrip() {
        let dir = std::env::temp_dir();
        let path = dir.join("test_sp3_product_writer.sp3");
        let scs = fake_spacecraft_states();
        Sp3ProductWriter::new("L01")
            .write_from_states(&path, &scs)
            .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with("PL01")).count(), 5);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn oem_product_writer_roundtrip() {
        let dir = std::env::temp_dir();
        let path = dir.join("test_oem_product_writer.oem");
        let scs = fake_spacecraft_states();
        OemProductWriter::new("2024-001A", "MYSAT")
            .write_from_states(&path, &scs)
            .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("OBJECT_ID            = 2024-001A"));
        let _ = std::fs::remove_file(&path);
    }
}
