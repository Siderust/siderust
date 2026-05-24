//! # CCSDS OPM reader/writer (CCSDS 502.0-B-3)
//!
//! ## Scientific scope
//!
//! Parses and emits Orbit Parameter Messages in Keyword-Value Notation (KVN).
//! Supported blocks: header, `META_START`/`META_STOP`, and state/vector
//! keywords (`EPOCH`, `X`/`Y`/`Z`, `X_DOT`/`Y_DOT`/`Z_DOT`).
//!
//! ## Technical scope
//!
//! The types in this module are plain CCSDS message containers, not the typed
//! scientific orbit API. They preserve the field names and units used by the
//! wire format so readers/writers can round-trip OPM payloads losslessly.
//!
//! ## Wire-format note
//!
//! [`CartesianState`] and [`KeplerianElements`] directly mirror the CCSDS OPM
//! fields and therefore use raw km, km/s, and degree scalars. Convert into the
//! typed `qtty`/`tempoch`/`affn` model after parsing when you need scientific
//! computation.
//!
//! ## References
//!
//! - CCSDS 502.0-B-3: Orbit Data Messages, Blue Book (2019).

use super::FormatError;
use std::io::{BufRead, BufReader, Read, Write};

/// OPM state vector (Cartesian).
///
/// This struct directly mirrors the CCSDS OPM `X`/`Y`/`Z` and
/// `X_DOT`/`Y_DOT`/`Z_DOT` fields, so its coordinates remain raw km and km/s
/// wire-format scalars.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::CartesianState;
/// let s = CartesianState {
///     epoch: "2024-001T12:00:00.000".to_string(),
///     position_km: [7000.0, 0.0, 0.0],
///     velocity_km_s: [0.0, 7.5, 0.0],
/// };
/// assert_eq!(s.position_km[0], 7000.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CartesianState {
    /// Epoch in ISO 8601 format.
    pub epoch: String,
    /// Position X, Y, Z in km.
    pub position_km: [f64; 3],
    /// Velocity Vx, Vy, Vz in km/s.
    pub velocity_km_s: [f64; 3],
}

/// OPM Keplerian elements section.
///
/// This struct directly mirrors the CCSDS OPM Keplerian-element fields, so the
/// values remain raw km and degree wire-format scalars.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::KeplerianElements;
/// let k = KeplerianElements {
///     semi_major_axis_km: 7000.0,
///     eccentricity: 0.001,
///     inclination_deg: 51.6,
///     ra_of_asc_node_deg: 247.0,
///     arg_of_pericenter_deg: 130.0,
///     true_anomaly_deg: 325.0,
/// };
/// assert!(k.semi_major_axis_km > 0.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct KeplerianElements {
    /// Semi-major axis in km.
    pub semi_major_axis_km: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Inclination in degrees.
    pub inclination_deg: f64,
    /// RAAN in degrees.
    pub ra_of_asc_node_deg: f64,
    /// Argument of perigee in degrees.
    pub arg_of_pericenter_deg: f64,
    /// True anomaly in degrees.
    pub true_anomaly_deg: f64,
}

/// OPM metadata block.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::OpmMetadata;
/// let m = OpmMetadata {
///     object_name: "TEST-SAT".to_string(),
///     object_id: "2024-001A".to_string(),
///     center_name: "EARTH".to_string(),
///     ref_frame: "EME2000".to_string(),
///     time_system: "UTC".to_string(),
/// };
/// assert_eq!(m.object_name, "TEST-SAT");
/// ```
#[derive(Debug, Clone)]
pub struct OpmMetadata {
    /// Object name.
    pub object_name: String,
    /// Object ID (COSPAR or user).
    pub object_id: String,
    /// Center body.
    pub center_name: String,
    /// Reference frame.
    pub ref_frame: String,
    /// Time system.
    pub time_system: String,
}

/// Parsed CCSDS OPM message.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::{OpmMessage, OpmMetadata, CartesianState};
/// let msg = OpmMessage {
///     metadata: OpmMetadata {
///         object_name: "SAT".to_string(),
///         object_id: "2024-001A".to_string(),
///         center_name: "EARTH".to_string(),
///         ref_frame: "EME2000".to_string(),
///         time_system: "UTC".to_string(),
///     },
///     state: CartesianState {
///         epoch: "2024-001T00:00:00.000".to_string(),
///         position_km: [7000.0, 0.0, 0.0],
///         velocity_km_s: [0.0, 7.5, 0.0],
///     },
///     keplerian: None,
/// };
/// assert_eq!(msg.state.position_km[0], 7000.0);
/// ```
#[derive(Debug, Clone)]
pub struct OpmMessage {
    /// Metadata block.
    pub metadata: OpmMetadata,
    /// Cartesian state (mandatory).
    pub state: CartesianState,
    /// Keplerian elements (optional).
    pub keplerian: Option<KeplerianElements>,
}

/// Read a CCSDS OPM KVN message from a byte source.
///
/// Parses `CCSDS_OPM_VERS`, `META_START`/`META_STOP`, state vector keywords,
/// and optionally the `KEPLERIAN` block.
///
/// # Errors
///
/// Returns [`FormatError::Format`] if mandatory fields are absent or malformed.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::read_opm;
///
/// let data = b"CCSDS_OPM_VERS = 2.0\nMETA_START\nOBJECT_NAME = SAT\n\
///              OBJECT_ID = 2024-001A\nCENTER_NAME = EARTH\nREF_FRAME = EME2000\n\
///              TIME_SYSTEM = UTC\nMETA_STOP\n\
///              EPOCH = 2024-001T12:00:00.000\nX = 7000.0\nY = 0.0\nZ = 0.0\n\
///              X_DOT = 0.0\nY_DOT = 7.5\nZ_DOT = 0.0\n";
/// let msg = read_opm(&data[..]).unwrap();
/// assert_eq!(msg.metadata.object_name, "SAT");
/// ```
pub fn read_opm<R: Read>(reader: R) -> Result<OpmMessage, FormatError> {
    let buf = BufReader::new(reader);
    let mut kv: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    for result in buf.lines() {
        let line = result.map_err(FormatError::Io)?;
        let line = line.trim().to_string();
        if line.is_empty()
            || line.starts_with("COMMENT")
            || line == "META_START"
            || line == "META_STOP"
        {
            continue;
        }
        if let Some(eq) = line.find('=') {
            let key = line[..eq].trim().to_string();
            let val = line[eq + 1..].trim().to_string();
            kv.insert(key, val);
        }
    }

    let get = |k: &str| -> Result<String, FormatError> {
        kv.get(k)
            .cloned()
            .ok_or_else(|| FormatError::Format(format!("OPM: missing field {k}")))
    };
    let getf = |k: &str| -> Result<f64, FormatError> {
        let v = kv
            .get(k)
            .ok_or_else(|| FormatError::Format(format!("OPM: missing field {k}")))?;
        v.parse::<f64>()
            .map_err(|_| FormatError::Format(format!("OPM: cannot parse {k}={v:?}")))
    };

    let metadata = OpmMetadata {
        object_name: get("OBJECT_NAME")?,
        object_id: get("OBJECT_ID").unwrap_or_default(),
        center_name: get("CENTER_NAME").unwrap_or_else(|_| "EARTH".to_string()),
        ref_frame: get("REF_FRAME").unwrap_or_else(|_| "EME2000".to_string()),
        time_system: get("TIME_SYSTEM").unwrap_or_else(|_| "UTC".to_string()),
    };

    let epoch = get("EPOCH")?;
    let state = CartesianState {
        epoch,
        position_km: [getf("X")?, getf("Y")?, getf("Z")?],
        velocity_km_s: [getf("X_DOT")?, getf("Y_DOT")?, getf("Z_DOT")?],
    };

    let keplerian = if kv.contains_key("SEMI_MAJOR_AXIS") {
        Some(KeplerianElements {
            semi_major_axis_km: getf("SEMI_MAJOR_AXIS")?,
            eccentricity: getf("ECCENTRICITY")?,
            inclination_deg: getf("INCLINATION")?,
            ra_of_asc_node_deg: getf("RA_OF_ASC_NODE")?,
            arg_of_pericenter_deg: getf("ARG_OF_PERICENTER")?,
            true_anomaly_deg: getf("TRUE_ANOMALY")?,
        })
    } else {
        None
    };

    Ok(OpmMessage {
        metadata,
        state,
        keplerian,
    })
}

/// Write a CCSDS OPM message in KVN format.
///
/// Emits the header, `META_START`/`META_STOP` block, and state vector keywords.
///
/// # Errors
///
/// Returns [`FormatError::Io`] on write failure.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::opm::{write_opm, read_opm, OpmMessage, OpmMetadata, CartesianState};
///
/// let msg = OpmMessage {
///     metadata: OpmMetadata {
///         object_name: "SAT".to_string(),
///         object_id: "2024-001A".to_string(),
///         center_name: "EARTH".to_string(),
///         ref_frame: "EME2000".to_string(),
///         time_system: "UTC".to_string(),
///     },
///     state: CartesianState {
///         epoch: "2024-001T12:00:00.000".to_string(),
///         position_km: [7000.0, 0.0, 0.0],
///         velocity_km_s: [0.0, 7.5, 0.0],
///     },
///     keplerian: None,
/// };
/// let mut buf = Vec::new();
/// write_opm(&mut buf, &msg).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.contains("OBJECT_NAME"));
/// ```
pub fn write_opm<W: Write>(w: &mut W, msg: &OpmMessage) -> Result<(), FormatError> {
    writeln!(w, "CCSDS_OPM_VERS = 2.0").map_err(FormatError::Io)?;
    writeln!(w, "META_START").map_err(FormatError::Io)?;
    writeln!(w, "OBJECT_NAME          = {}", msg.metadata.object_name).map_err(FormatError::Io)?;
    writeln!(w, "OBJECT_ID            = {}", msg.metadata.object_id).map_err(FormatError::Io)?;
    writeln!(w, "CENTER_NAME          = {}", msg.metadata.center_name).map_err(FormatError::Io)?;
    writeln!(w, "REF_FRAME            = {}", msg.metadata.ref_frame).map_err(FormatError::Io)?;
    writeln!(w, "TIME_SYSTEM          = {}", msg.metadata.time_system).map_err(FormatError::Io)?;
    writeln!(w, "META_STOP").map_err(FormatError::Io)?;
    writeln!(w).map_err(FormatError::Io)?;
    writeln!(w, "EPOCH                = {}", msg.state.epoch).map_err(FormatError::Io)?;
    writeln!(w, "X                    = {:.6}", msg.state.position_km[0])
        .map_err(FormatError::Io)?;
    writeln!(w, "Y                    = {:.6}", msg.state.position_km[1])
        .map_err(FormatError::Io)?;
    writeln!(w, "Z                    = {:.6}", msg.state.position_km[2])
        .map_err(FormatError::Io)?;
    writeln!(
        w,
        "X_DOT                = {:.6}",
        msg.state.velocity_km_s[0]
    )
    .map_err(FormatError::Io)?;
    writeln!(
        w,
        "Y_DOT                = {:.6}",
        msg.state.velocity_km_s[1]
    )
    .map_err(FormatError::Io)?;
    writeln!(
        w,
        "Z_DOT                = {:.6}",
        msg.state.velocity_km_s[2]
    )
    .map_err(FormatError::Io)?;
    if let Some(k) = &msg.keplerian {
        writeln!(w).map_err(FormatError::Io)?;
        writeln!(w, "SEMI_MAJOR_AXIS      = {:.6}", k.semi_major_axis_km)
            .map_err(FormatError::Io)?;
        writeln!(w, "ECCENTRICITY         = {:.7}", k.eccentricity).map_err(FormatError::Io)?;
        writeln!(w, "INCLINATION          = {:.4}", k.inclination_deg).map_err(FormatError::Io)?;
        writeln!(w, "RA_OF_ASC_NODE       = {:.4}", k.ra_of_asc_node_deg)
            .map_err(FormatError::Io)?;
        writeln!(w, "ARG_OF_PERICENTER    = {:.4}", k.arg_of_pericenter_deg)
            .map_err(FormatError::Io)?;
        writeln!(w, "TRUE_ANOMALY         = {:.4}", k.true_anomaly_deg).map_err(FormatError::Io)?;
    }
    Ok(())
}
