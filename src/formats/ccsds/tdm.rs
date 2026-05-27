//! # CCSDS TDM reader/writer (CCSDS 503.0-B-2)
//!
//! Parses and emits Tracking Data Messages in Keyword-Value Notation (KVN).
//! Supported blocks: header, `META_START`/`META_STOP`,
//! `DATA_START`/`DATA_STOP`.
//!
//! ## References
//!
//! - CCSDS 503.0-B-2: Tracking Data Message, Blue Book (2007).

use super::FormatError;
use std::io::{BufRead, BufReader, Read, Write};

/// TDM observation type.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::ObservationType;
/// assert_eq!(ObservationType::Range, ObservationType::Range);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum ObservationType {
    /// Two-way range measurement.
    Range,
    /// Doppler frequency shift (instantaneous).
    Doppler,
    /// Angle 1 (azimuth or right ascension).
    Angle1,
    /// Angle 2 (elevation or declination).
    Angle2,
    /// Other observable keyword.
    Other(String),
}

/// A single TDM observation record.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::{ObservationData, ObservationType};
/// let obs = ObservationData {
///     obs_type: ObservationType::Range,
///     epoch: "2024-001T12:00:00".to_string(),
///     value: 23000.0,
/// };
/// assert_eq!(obs.value, 23000.0);
/// ```
#[derive(Debug, Clone)]
pub struct ObservationData {
    /// Observation type.
    pub obs_type: ObservationType,
    /// Epoch in ISO 8601 format.
    pub epoch: String,
    /// Observation value (units depend on type).
    pub value: f64,
}

/// TDM metadata block.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::TdmMetadata;
/// let m = TdmMetadata {
///     participants: vec!["STATION_A".to_string()],
///     mode: "SEQUENTIAL".to_string(),
///     path: "1,2,1".to_string(),
///     time_system: "UTC".to_string(),
/// };
/// assert_eq!(m.mode, "SEQUENTIAL");
/// ```
#[derive(Debug, Clone)]
pub struct TdmMetadata {
    /// Participant identifiers.
    pub participants: Vec<String>,
    /// Mode (e.g. `"SEQUENTIAL"`).
    pub mode: String,
    /// Path string (e.g. `"1,2,1"`).
    pub path: String,
    /// Time system.
    pub time_system: String,
}

/// Parsed CCSDS TDM message.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::TdmMessage;
/// let msg = TdmMessage::default();
/// assert!(msg.observations.is_empty());
/// ```
#[derive(Debug, Default)]
pub struct TdmMessage {
    /// Metadata block(s).
    pub metadata: Vec<TdmMetadata>,
    /// Observation records.
    pub observations: Vec<ObservationData>,
}

/// Read a CCSDS TDM KVN message from a byte source.
///
/// Parses `CCSDS_TDM_VERS`, `META_START`/`META_STOP`,
/// and `DATA_START`/`DATA_STOP` blocks.
///
/// # Errors
///
/// Returns [`FormatError::Format`] for missing or malformed fields.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::read_tdm;
///
/// let data = b"CCSDS_TDM_VERS = 1.0\nMETA_START\nPARTICIPANT_1 = STA\n\
///              MODE = SEQUENTIAL\nPATH = 1,2,1\nTIME_SYSTEM = UTC\nMETA_STOP\n\
///              DATA_START\nRANGE = 2024-001T12:00:00 : 23000.0\nDATA_STOP\n";
/// let msg = read_tdm(&data[..]).unwrap();
/// assert_eq!(msg.observations.len(), 1);
/// ```
pub fn read_tdm<R: Read>(reader: R) -> Result<TdmMessage, FormatError> {
    let buf = BufReader::new(reader);
    let mut msg = TdmMessage::default();

    enum State {
        None,
        Meta,
        Data,
    }
    let mut state = State::None;
    let mut current_meta = new_meta();

    for result in buf.lines() {
        let line = result.map_err(FormatError::Io)?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with("COMMENT") || line.starts_with("CCSDS_TDM") {
            continue;
        }
        match line.as_str() {
            "META_START" => {
                state = State::Meta;
                current_meta = new_meta();
                continue;
            }
            "META_STOP" => {
                msg.metadata.push(current_meta.clone());
                current_meta = new_meta();
                state = State::None;
                continue;
            }
            "DATA_START" => {
                state = State::Data;
                continue;
            }
            "DATA_STOP" => {
                state = State::None;
                continue;
            }
            _ => {}
        }

        if let Some(eq) = line.find('=') {
            let key = line[..eq].trim().to_string();
            let val = line[eq + 1..].trim().to_string();
            match state {
                State::Meta => {
                    if key.starts_with("PARTICIPANT_") {
                        current_meta.participants.push(val);
                    } else if key == "MODE" {
                        current_meta.mode = val;
                    } else if key == "PATH" {
                        current_meta.path = val;
                    } else if key == "TIME_SYSTEM" {
                        current_meta.time_system = val;
                    }
                }
                State::Data => {
                    // val format: "epoch : value" — split on last " : " to avoid
                    // colons in ISO 8601 epoch strings (e.g. "2024-001T12:00:00 : 23000.0")
                    let sep = " : ";
                    let colon = val.rfind(sep).ok_or_else(|| {
                        FormatError::Format(format!("TDM: data line missing ' : ': {key} = {val}"))
                    })?;
                    let epoch = val[..colon].trim().to_string();
                    let value_str = val[colon + sep.len()..].trim();
                    let value = value_str.parse::<f64>().map_err(|_| {
                        FormatError::Format(format!("TDM: cannot parse value: {value_str:?}"))
                    })?;
                    let obs_type = keyword_to_obs_type(&key);
                    msg.observations.push(ObservationData {
                        obs_type,
                        epoch,
                        value,
                    });
                }
                State::None => {}
            }
        }
    }

    Ok(msg)
}

fn new_meta() -> TdmMetadata {
    TdmMetadata {
        participants: Vec::new(),
        mode: String::new(),
        path: String::new(),
        time_system: String::new(),
    }
}

fn keyword_to_obs_type(kw: &str) -> ObservationType {
    match kw {
        "RANGE" => ObservationType::Range,
        "DOPPLER_INSTANTANEOUS" | "DOPPLER_INTEGRATED" => ObservationType::Doppler,
        "ANGLE_1" => ObservationType::Angle1,
        "ANGLE_2" => ObservationType::Angle2,
        other => ObservationType::Other(other.to_string()),
    }
}

fn obs_type_to_keyword(t: &ObservationType) -> &str {
    match t {
        ObservationType::Range => "RANGE",
        ObservationType::Doppler => "DOPPLER_INSTANTANEOUS",
        ObservationType::Angle1 => "ANGLE_1",
        ObservationType::Angle2 => "ANGLE_2",
        ObservationType::Other(s) => s.as_str(),
    }
}

/// Write a CCSDS TDM message in KVN format.
///
/// Emits the header, `META_START`/`META_STOP` block(s), and
/// `DATA_START`/`DATA_STOP` block.
///
/// # Errors
///
/// Returns [`FormatError::Io`] on write failure.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::tdm::{write_tdm, TdmMessage, TdmMetadata, ObservationData, ObservationType};
///
/// let mut msg = TdmMessage::default();
/// msg.metadata.push(TdmMetadata {
///     participants: vec!["STA".to_string()],
///     mode: "SEQUENTIAL".to_string(),
///     path: "1,2,1".to_string(),
///     time_system: "UTC".to_string(),
/// });
/// msg.observations.push(ObservationData {
///     obs_type: ObservationType::Range,
///     epoch: "2024-001T12:00:00".to_string(),
///     value: 23000.0,
/// });
/// let mut buf = Vec::new();
/// write_tdm(&mut buf, &msg).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.contains("RANGE"));
/// ```
pub fn write_tdm<W: Write>(w: &mut W, msg: &TdmMessage) -> Result<(), FormatError> {
    writeln!(w, "CCSDS_TDM_VERS = 1.0").map_err(FormatError::Io)?;
    for (i, meta) in msg.metadata.iter().enumerate() {
        if i > 0 {
            writeln!(w).map_err(FormatError::Io)?;
        }
        writeln!(w, "META_START").map_err(FormatError::Io)?;
        for (j, p) in meta.participants.iter().enumerate() {
            writeln!(w, "PARTICIPANT_{} = {}", j + 1, p).map_err(FormatError::Io)?;
        }
        if !meta.mode.is_empty() {
            writeln!(w, "MODE             = {}", meta.mode).map_err(FormatError::Io)?;
        }
        if !meta.path.is_empty() {
            writeln!(w, "PATH             = {}", meta.path).map_err(FormatError::Io)?;
        }
        if !meta.time_system.is_empty() {
            writeln!(w, "TIME_SYSTEM      = {}", meta.time_system).map_err(FormatError::Io)?;
        }
        writeln!(w, "META_STOP").map_err(FormatError::Io)?;
    }
    writeln!(w).map_err(FormatError::Io)?;
    writeln!(w, "DATA_START").map_err(FormatError::Io)?;
    for obs in &msg.observations {
        writeln!(
            w,
            "{} = {} : {}",
            obs_type_to_keyword(&obs.obs_type),
            obs.epoch,
            obs.value
        )
        .map_err(FormatError::Io)?;
    }
    writeln!(w, "DATA_STOP").map_err(FormatError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_tdm_text() -> &'static str {
        "CCSDS_TDM_VERS = 1.0\n\
         META_START\n\
         PARTICIPANT_1 = STATION_A\n\
         PARTICIPANT_2 = SAT1\n\
         MODE = SEQUENTIAL\n\
         PATH = 1,2,1\n\
         TIME_SYSTEM = UTC\n\
         META_STOP\n\
         DATA_START\n\
         RANGE = 2024-001T12:00:00 : 23000.0\n\
         DATA_STOP\n"
    }

    #[test]
    fn read_tdm_basic_range() {
        let msg = read_tdm(base_tdm_text().as_bytes()).unwrap();
        assert_eq!(msg.observations.len(), 1);
        assert_eq!(msg.observations[0].obs_type, ObservationType::Range);
        assert_eq!(msg.observations[0].value, 23000.0);
        assert_eq!(msg.observations[0].epoch, "2024-001T12:00:00");
    }

    #[test]
    fn read_tdm_all_obs_types() {
        let text = "CCSDS_TDM_VERS = 1.0\n\
                    META_START\nPARTICIPANT_1 = S\nMETA_STOP\n\
                    DATA_START\n\
                    RANGE = 2024-001T00:00:00 : 1.0\n\
                    DOPPLER_INSTANTANEOUS = 2024-001T00:01:00 : 2.0\n\
                    ANGLE_1 = 2024-001T00:02:00 : 3.0\n\
                    ANGLE_2 = 2024-001T00:03:00 : 4.0\n\
                    MY_CUSTOM_OBS = 2024-001T00:04:00 : 5.0\n\
                    DATA_STOP\n";
        let msg = read_tdm(text.as_bytes()).unwrap();
        assert_eq!(msg.observations.len(), 5);
        assert_eq!(msg.observations[1].obs_type, ObservationType::Doppler);
        assert_eq!(msg.observations[2].obs_type, ObservationType::Angle1);
        assert_eq!(msg.observations[3].obs_type, ObservationType::Angle2);
        assert!(
            matches!(&msg.observations[4].obs_type, ObservationType::Other(s) if s == "MY_CUSTOM_OBS")
        );
    }

    #[test]
    fn read_tdm_multiple_metadata_blocks() {
        let text = "CCSDS_TDM_VERS = 1.0\n\
                    META_START\nPARTICIPANT_1 = A\nMETA_STOP\n\
                    META_START\nPARTICIPANT_1 = B\nMETA_STOP\n\
                    DATA_START\nRANGE = 2024-001T00:00:00 : 10.0\nDATA_STOP\n";
        let msg = read_tdm(text.as_bytes()).unwrap();
        assert_eq!(msg.metadata.len(), 2);
        assert_eq!(msg.metadata[0].participants[0], "A");
        assert_eq!(msg.metadata[1].participants[0], "B");
    }

    #[test]
    fn read_tdm_missing_separator_is_error() {
        let text = "CCSDS_TDM_VERS = 1.0\n\
                    META_START\nPARTICIPANT_1 = S\nMETA_STOP\n\
                    DATA_START\n\
                    RANGE = 2024-001T00:00:00_NO_SEPARATOR\n\
                    DATA_STOP\n";
        assert!(read_tdm(text.as_bytes()).is_err());
    }

    #[test]
    fn read_tdm_invalid_float_is_error() {
        let text = "CCSDS_TDM_VERS = 1.0\n\
                    META_START\nPARTICIPANT_1 = S\nMETA_STOP\n\
                    DATA_START\n\
                    RANGE = 2024-001T00:00:00 : not_a_float\n\
                    DATA_STOP\n";
        assert!(read_tdm(text.as_bytes()).is_err());
    }

    #[test]
    fn write_tdm_roundtrip() {
        let mut msg = TdmMessage::default();
        msg.metadata.push(TdmMetadata {
            participants: vec!["STATION_A".to_string(), "SAT1".to_string()],
            mode: "SEQUENTIAL".to_string(),
            path: "1,2,1".to_string(),
            time_system: "UTC".to_string(),
        });
        msg.observations.push(ObservationData {
            obs_type: ObservationType::Range,
            epoch: "2024-001T12:00:00".to_string(),
            value: 23000.0,
        });
        msg.observations.push(ObservationData {
            obs_type: ObservationType::Angle1,
            epoch: "2024-001T12:01:00".to_string(),
            value: 45.5,
        });
        let mut buf = Vec::new();
        write_tdm(&mut buf, &msg).unwrap();
        let text = String::from_utf8(buf.clone()).unwrap();
        assert!(text.contains("RANGE"));
        assert!(text.contains("ANGLE_1"));
        let parsed = read_tdm(buf.as_slice()).unwrap();
        assert_eq!(parsed.observations.len(), 2);
        assert_eq!(parsed.observations[0].value, 23000.0);
    }
}
