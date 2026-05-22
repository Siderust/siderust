//! # CCSDS OEM writer
//!
//! ## Scientific scope
//!
//! The Orbit Ephemeris Message is a standard interchange product for
//! spacecraft state histories. This module writes the subset needed by the
//! current POD pipeline: a deterministic ASCII KVN ephemeris segment
//! without covariance.
//!
//! Its scientific interpretation comes from the supplied state history and
//! metadata. No interpolation, force modelling, or covariance transport is
//! introduced during writing.
//!
//! ## Technical scope
//!
//! The public items are `OemMetadata` and `write_oem`. Callers provide a
//! stream of `OemState` values already expressed in the desired frame and
//! centre, and the writer serializes them into the CCSDS text layout.
//!
//! Reading OEM files and handling richer message features such as maneuver
//! blocks or covariance segments are out of scope here.
//!
//! ## References
//!
//! - Consultative Committee for Space Data Systems. (2010). Orbit Data
//!   Messages, CCSDS 502.0-B-2 / 502.0-B-3.
//! - Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications
//!   (4th ed.). Microcosm Press.
use super::FormatError;
use std::io::Write;

/// A spacecraft state as stored in a CCSDS OEM file.
///
/// The reference frame and time system are metadata-level: they are recorded
/// in the containing [`OemMetadata`] rather than in the type, because OEM
/// supports multiple frames (EME2000, GCRF, ITRF, …) chosen at write time.
/// All numeric fields use the conventional OEM units: km for position,
/// km/s for velocity, Julian Date (any time system) for epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OemState {
    /// Epoch expressed as a Julian Date in the time system given by the metadata.
    pub epoch_jd: f64,
    /// X, Y, Z position components (km).
    pub position_km: [f64; 3],
    /// VX, VY, VZ velocity components (km/s).
    pub velocity_km_s: [f64; 3],
}

impl OemState {
    /// Construct an `OemState` from a Julian Date and component arrays.
    pub fn new(epoch_jd: f64, position_km: [f64; 3], velocity_km_s: [f64; 3]) -> Self {
        Self {
            epoch_jd,
            position_km,
            velocity_km_s,
        }
    }
}

/// Metadata for a single OEM segment.
#[derive(Debug, Clone)]
pub struct OemMetadata {
    /// Object identifier (free-form).
    pub object_id: String,
    /// Object name (free-form).
    pub object_name: String,
    /// Inertial reference frame string, e.g. "EME2000" or "GCRF".
    pub ref_frame: String,
    /// Time system string, e.g. "TT" or "UTC".
    pub time_system: String,
    /// Centre body, e.g. "EARTH".
    pub center_name: String,
}

/// Write a CCSDS OEM file from a sequence of states.
pub fn write_oem<W: Write>(
    w: &mut W,
    meta: &OemMetadata,
    states: &[OemState],
) -> Result<(), FormatError> {
    if states.is_empty() {
        return Err(FormatError::Format(
            "write_oem: no states to write".to_string(),
        ));
    }
    writeln!(w, "CCSDS_OEM_VERS = 3.0")?;
    writeln!(w, "CREATION_DATE  = 2026-01-01T00:00:00")?;
    writeln!(w, "ORIGINATOR     = SIDERUST-POD")?;
    writeln!(w)?;
    writeln!(w, "META_START")?;
    writeln!(w, "OBJECT_NAME          = {}", meta.object_name)?;
    writeln!(w, "OBJECT_ID            = {}", meta.object_id)?;
    writeln!(w, "CENTER_NAME          = {}", meta.center_name)?;
    writeln!(w, "REF_FRAME            = {}", meta.ref_frame)?;
    writeln!(w, "TIME_SYSTEM          = {}", meta.time_system)?;
    writeln!(
        w,
        "START_TIME           = {}",
        jd_to_iso8601(states[0].epoch_jd)
    )?;
    writeln!(
        w,
        "STOP_TIME            = {}",
        jd_to_iso8601(states.last().unwrap().epoch_jd)
    )?;
    writeln!(w, "META_STOP")?;
    writeln!(w)?;
    for s in states {
        writeln!(
            w,
            "{} {:.6} {:.6} {:.6} {:.9} {:.9} {:.9}",
            jd_to_iso8601(s.epoch_jd),
            s.position_km[0],
            s.position_km[1],
            s.position_km[2],
            s.velocity_km_s[0],
            s.velocity_km_s[1],
            s.velocity_km_s[2],
        )?;
    }
    Ok(())
}

/// Convert a Julian Date (any time scale) to an ISO-8601-ish string.
///
/// Conversion is "calendar date in proleptic Gregorian, fractional seconds";
/// no leap-second handling — the caller is responsible for the time system.
fn jd_to_iso8601(jd: f64) -> String {
    let jd_int = (jd + 0.5).floor() as i64;
    let frac = jd + 0.5 - jd_int as f64;

    let a = jd_int + 32_044;
    let b = (4 * a + 3) / 146_097;
    let c = a - (146_097 * b) / 4;
    let d = (4 * c + 3) / 1_461;
    let e = c - (1_461 * d) / 4;
    let m = (5 * e + 2) / 153;
    let day = (e - (153 * m + 2) / 5 + 1) as u32;
    let month = (m + 3 - 12 * (m / 10)) as u32;
    let year = (100 * b + d - 4_800 + (m / 10)) as i32;

    let total_seconds = frac * 86_400.0;
    let hours = (total_seconds / 3600.0).floor() as u32 % 24;
    let mins = ((total_seconds % 3600.0) / 60.0).floor() as u32;
    let secs = total_seconds % 60.0;

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:09.6}",
        year, month, day, hours, mins, secs
    )
}

/// A parsed OEM segment.
#[derive(Debug, Clone)]
pub struct OemSegment {
    /// Metadata block.
    pub metadata: OemMetadata,
    /// State history in the order they appear in the file.
    pub states: Vec<OemState>,
}

/// A parsed OEM file (header + one or more segments).
#[derive(Debug, Clone)]
pub struct OemFile {
    /// `CCSDS_OEM_VERS` value (e.g. `"3.0"`).
    pub version: String,
    /// `CREATION_DATE` header field, verbatim.
    pub creation_date: String,
    /// `ORIGINATOR` header field, verbatim.
    pub originator: String,
    /// One segment per `META_START`/`META_STOP` block in input order.
    pub segments: Vec<OemSegment>,
}

/// Read a CCSDS OEM (KVN, version 2.x/3.x) file from `r`.
///
/// The reader accepts the subset that [`write_oem`] produces *and* the
/// canonical CCSDS layout: arbitrary whitespace inside data lines, optional
/// `COMMENT` lines anywhere outside a data block, and one or more segments.
/// Covariance and maneuver sub-blocks are tolerated but skipped.
///
/// # Errors
///
/// Returns [`FormatError::Format`] for any malformed line, missing required
/// metadata field, unparseable epoch, or non-finite numeric component.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::oem::{read_oem, write_oem, OemMetadata};
/// use siderust::formats::ccsds::oem::OemState;
///
/// let states = vec![OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0])];
/// let meta = OemMetadata {
///     object_id: "1900-001A".into(),
///     object_name: "TEST".into(),
///     ref_frame: "EME2000".into(),
///     time_system: "TT".into(),
///     center_name: "EARTH".into(),
/// };
/// let mut buf = Vec::new();
/// write_oem(&mut buf, &meta, &states).unwrap();
/// let parsed = read_oem(&buf[..]).unwrap();
/// assert_eq!(parsed.segments.len(), 1);
/// assert_eq!(parsed.segments[0].states.len(), 1);
/// ```
pub fn read_oem<R: std::io::Read>(mut r: R) -> Result<OemFile, FormatError> {
    let mut text = String::new();
    r.read_to_string(&mut text)
        .map_err(|e| FormatError::Format(format!("read_oem: io error: {e}")))?;

    let mut version: Option<String> = None;
    let mut creation_date: Option<String> = None;
    let mut originator: Option<String> = None;
    let mut segments: Vec<OemSegment> = Vec::new();

    enum State {
        Header,
        Meta(MetaBuilder),
        Data(OemSegment),
        SkipBlock,
    }

    #[derive(Default)]
    struct MetaBuilder {
        object_id: Option<String>,
        object_name: Option<String>,
        ref_frame: Option<String>,
        time_system: Option<String>,
        center_name: Option<String>,
    }

    let mut state = State::Header;
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with("COMMENT") {
            continue;
        }
        match &mut state {
            State::Header => {
                if line == "META_START" {
                    state = State::Meta(MetaBuilder::default());
                    continue;
                }
                let (k, v) = split_kv(line, lineno)?;
                match k {
                    "CCSDS_OEM_VERS" => version = Some(v.into()),
                    "CREATION_DATE" => creation_date = Some(v.into()),
                    "ORIGINATOR" => originator = Some(v.into()),
                    _ => {}
                }
            }
            State::Meta(mb) => {
                if line == "META_STOP" {
                    let MetaBuilder {
                        object_id,
                        object_name,
                        ref_frame,
                        time_system,
                        center_name,
                    } = std::mem::take(mb);
                    let need = |o: Option<String>, name: &str| -> Result<String, FormatError> {
                        o.ok_or_else(|| {
                            FormatError::Format(format!(
                                "read_oem: META block missing required field {name}"
                            ))
                        })
                    };
                    let metadata = OemMetadata {
                        object_id: need(object_id, "OBJECT_ID")?,
                        object_name: need(object_name, "OBJECT_NAME")?,
                        ref_frame: need(ref_frame, "REF_FRAME")?,
                        time_system: need(time_system, "TIME_SYSTEM")?,
                        center_name: need(center_name, "CENTER_NAME")?,
                    };
                    state = State::Data(OemSegment {
                        metadata,
                        states: Vec::new(),
                    });
                    continue;
                }
                let (k, v) = split_kv(line, lineno)?;
                match k {
                    "OBJECT_ID" => mb.object_id = Some(v.into()),
                    "OBJECT_NAME" => mb.object_name = Some(v.into()),
                    "REF_FRAME" => mb.ref_frame = Some(v.into()),
                    "TIME_SYSTEM" => mb.time_system = Some(v.into()),
                    "CENTER_NAME" => mb.center_name = Some(v.into()),
                    // START_TIME / STOP_TIME / USEABLE_*: informational, ignored.
                    _ => {}
                }
            }
            State::Data(seg) => {
                if line == "META_START" {
                    let done = std::mem::replace(
                        seg,
                        OemSegment {
                            metadata: seg.metadata.clone(),
                            states: Vec::new(),
                        },
                    );
                    segments.push(done);
                    state = State::Meta(MetaBuilder::default());
                    continue;
                }
                if line == "COVARIANCE_START" || line == "MAN_START" {
                    state = State::SkipBlock;
                    continue;
                }
                seg.states.push(parse_data_line(line, lineno)?);
            }
            State::SkipBlock => {
                if line == "COVARIANCE_STOP" || line == "MAN_STOP" {
                    // Re-enter the most recent Data segment.
                    if let Some(last) = segments.last_mut() {
                        let resumed = std::mem::replace(
                            last,
                            OemSegment {
                                metadata: last.metadata.clone(),
                                states: Vec::new(),
                            },
                        );
                        state = State::Data(resumed);
                    } else {
                        state = State::Header;
                    }
                }
            }
        }
    }
    if let State::Data(seg) = state {
        segments.push(seg);
    }
    if segments.is_empty() {
        return Err(FormatError::Format(
            "read_oem: file contains no segments".into(),
        ));
    }
    Ok(OemFile {
        version: version.unwrap_or_else(|| "unknown".into()),
        creation_date: creation_date.unwrap_or_default(),
        originator: originator.unwrap_or_default(),
        segments,
    })
}

fn split_kv(line: &str, lineno: usize) -> Result<(&str, &str), FormatError> {
    let (k, v) = line.split_once('=').ok_or_else(|| {
        FormatError::Format(format!(
            "read_oem: line {}: expected KEY = VALUE",
            lineno + 1
        ))
    })?;
    Ok((k.trim(), v.trim()))
}

fn parse_data_line(line: &str, lineno: usize) -> Result<OemState, FormatError> {
    let mut it = line.split_whitespace();
    let epoch_str = it.next().ok_or_else(|| {
        FormatError::Format(format!("read_oem: line {}: empty data row", lineno + 1))
    })?;
    let mut next_f64 = |what: &str| -> Result<f64, FormatError> {
        let tok = it.next().ok_or_else(|| {
            FormatError::Format(format!("read_oem: line {}: missing {what}", lineno + 1))
        })?;
        let v: f64 = tok.parse().map_err(|e| {
            FormatError::Format(format!(
                "read_oem: line {}: cannot parse {what} ({tok}): {e}",
                lineno + 1
            ))
        })?;
        if !v.is_finite() {
            return Err(FormatError::Format(format!(
                "read_oem: line {}: {what} not finite",
                lineno + 1
            )));
        }
        Ok(v)
    };
    let x = next_f64("X")?;
    let y = next_f64("Y")?;
    let z = next_f64("Z")?;
    let vx = next_f64("VX")?;
    let vy = next_f64("VY")?;
    let vz = next_f64("VZ")?;
    let jd = iso8601_to_jd(epoch_str).ok_or_else(|| {
        FormatError::Format(format!(
            "read_oem: line {}: cannot parse epoch {epoch_str}",
            lineno + 1
        ))
    })?;
    Ok(OemState::new(jd, [x, y, z], [vx, vy, vz]))
}

/// Inverse of [`jd_to_iso8601`]. Returns `None` on malformed input.
///
/// Accepts `YYYY-MM-DDTHH:MM:SS[.ffffff]`. The conversion is purely calendrical
/// (no leap-second handling); the caller owns the time-system semantics.
fn iso8601_to_jd(s: &str) -> Option<f64> {
    let (date, time) = s.split_once('T')?;
    let mut date_parts = date.split('-');
    let year: i32 = date_parts.next()?.parse().ok()?;
    let month: u32 = date_parts.next()?.parse().ok()?;
    let day: u32 = date_parts.next()?.parse().ok()?;
    let mut time_parts = time.split(':');
    let hour: u32 = time_parts.next()?.parse().ok()?;
    let minute: u32 = time_parts.next()?.parse().ok()?;
    let second: f64 = time_parts.next()?.parse().ok()?;

    // Proleptic Gregorian calendar → JDN at 0h UT.
    let m = month as i64;
    let y = year as i64;
    let a = (14 - m) / 12;
    let yy = y + 4800 - a;
    let mm = m + 12 * a - 3;
    let jdn = day as i64 + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32_045;
    let day_frac = (hour as f64 * 3600.0 + minute as f64 * 60.0 + second) / 86_400.0;
    Some(jdn as f64 - 0.5 + day_frac)
}

// ─── OEM XML (CCSDS 502.0-B-3 §6) ──────────────────────────────────────────

/// Write a CCSDS OEM XML (`<oem>` root) document equivalent to the KVN
/// produced by [`write_oem`].
///
/// One `<segment>` is emitted with `<metadata>` and `<data>/<stateVector>`
/// children. Epochs use the same ISO-8601 calendrical conversion as
/// [`write_oem`].
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::oem::{write_oem_xml, OemMetadata};
/// use siderust::formats::ccsds::oem::OemState;
///
/// let states = vec![OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0])];
/// let meta = OemMetadata {
///     object_id: "1900-001A".into(),
///     object_name: "TEST".into(),
///     ref_frame: "EME2000".into(),
///     time_system: "TT".into(),
///     center_name: "EARTH".into(),
/// };
/// let mut buf = Vec::new();
/// write_oem_xml(&mut buf, &meta, &states).unwrap();
/// let txt = std::str::from_utf8(&buf).unwrap();
/// assert!(txt.contains("<oem"));
/// assert!(txt.contains("EME2000"));
/// ```
pub fn write_oem_xml<W: std::io::Write>(
    w: &mut W,
    meta: &OemMetadata,
    states: &[OemState],
) -> Result<(), FormatError> {
    if states.is_empty() {
        return Err(FormatError::Format(
            "write_oem_xml: no states to write".into(),
        ));
    }
    let start = jd_to_iso8601(states[0].epoch_jd);
    let stop = jd_to_iso8601(states.last().unwrap().epoch_jd);
    writeln!(w, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")?;
    writeln!(w, "<oem id=\"CCSDS_OEM_VERS\" version=\"3.0\">")?;
    writeln!(w, "  <header>")?;
    writeln!(w, "    <CREATION_DATE>2026-01-01T00:00:00</CREATION_DATE>")?;
    writeln!(w, "    <ORIGINATOR>SIDERUST-POD</ORIGINATOR>")?;
    writeln!(w, "  </header>")?;
    writeln!(w, "  <body>")?;
    writeln!(w, "    <segment>")?;
    writeln!(w, "      <metadata>")?;
    writeln!(
        w,
        "        <OBJECT_NAME>{}</OBJECT_NAME>",
        xml_escape(&meta.object_name)
    )?;
    writeln!(
        w,
        "        <OBJECT_ID>{}</OBJECT_ID>",
        xml_escape(&meta.object_id)
    )?;
    writeln!(
        w,
        "        <CENTER_NAME>{}</CENTER_NAME>",
        xml_escape(&meta.center_name)
    )?;
    writeln!(
        w,
        "        <REF_FRAME>{}</REF_FRAME>",
        xml_escape(&meta.ref_frame)
    )?;
    writeln!(
        w,
        "        <TIME_SYSTEM>{}</TIME_SYSTEM>",
        xml_escape(&meta.time_system)
    )?;
    writeln!(w, "        <START_TIME>{}</START_TIME>", start)?;
    writeln!(w, "        <STOP_TIME>{}</STOP_TIME>", stop)?;
    writeln!(w, "      </metadata>")?;
    writeln!(w, "      <data>")?;
    for s in states {
        writeln!(w, "        <stateVector>")?;
        writeln!(w, "          <EPOCH>{}</EPOCH>", jd_to_iso8601(s.epoch_jd))?;
        writeln!(w, "          <X units=\"km\">{:.6}</X>", s.position_km[0])?;
        writeln!(w, "          <Y units=\"km\">{:.6}</Y>", s.position_km[1])?;
        writeln!(w, "          <Z units=\"km\">{:.6}</Z>", s.position_km[2])?;
        writeln!(
            w,
            "          <X_DOT units=\"km/s\">{:.9}</X_DOT>",
            s.velocity_km_s[0]
        )?;
        writeln!(
            w,
            "          <Y_DOT units=\"km/s\">{:.9}</Y_DOT>",
            s.velocity_km_s[1]
        )?;
        writeln!(
            w,
            "          <Z_DOT units=\"km/s\">{:.9}</Z_DOT>",
            s.velocity_km_s[2]
        )?;
        writeln!(w, "        </stateVector>")?;
    }
    writeln!(w, "      </data>")?;
    writeln!(w, "    </segment>")?;
    writeln!(w, "  </body>")?;
    writeln!(w, "</oem>")?;
    Ok(())
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Read a CCSDS OEM XML (`<oem>` root) document.
///
/// The reader is the symmetric inverse of [`write_oem_xml`] and additionally
/// tolerates whitespace, comments, and namespace prefixes. It uses
/// `quick_xml` as the underlying parser.
///
/// # Errors
///
/// Returns [`FormatError::Format`] for any malformed XML, missing required
/// metadata field, unparseable epoch, or non-finite numeric component.
///
/// # Examples
///
/// ```
/// use siderust::formats::ccsds::oem::{read_oem_xml, write_oem_xml, OemMetadata};
/// use siderust::formats::ccsds::oem::OemState;
///
/// let states = vec![OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0])];
/// let meta = OemMetadata {
///     object_id: "1900-001A".into(),
///     object_name: "TEST".into(),
///     ref_frame: "EME2000".into(),
///     time_system: "TT".into(),
///     center_name: "EARTH".into(),
/// };
/// let mut buf = Vec::new();
/// write_oem_xml(&mut buf, &meta, &states).unwrap();
/// let parsed = read_oem_xml(&buf[..]).unwrap();
/// assert_eq!(parsed.segments.len(), 1);
/// assert_eq!(parsed.segments[0].states.len(), 1);
/// ```
pub fn read_oem_xml<R: std::io::Read>(mut r: R) -> Result<OemFile, FormatError> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let mut text = String::new();
    r.read_to_string(&mut text)
        .map_err(|e| FormatError::Format(format!("read_oem_xml: io: {e}")))?;
    let mut rdr = Reader::from_str(&text);
    rdr.config_mut().trim_text(true);

    #[derive(Default)]
    struct StateBuilder {
        epoch: Option<f64>,
        x: Option<f64>,
        y: Option<f64>,
        z: Option<f64>,
        vx: Option<f64>,
        vy: Option<f64>,
        vz: Option<f64>,
    }
    #[derive(Default)]
    struct SegBuilder {
        object_id: Option<String>,
        object_name: Option<String>,
        ref_frame: Option<String>,
        time_system: Option<String>,
        center_name: Option<String>,
        states: Vec<OemState>,
        cur: Option<StateBuilder>,
    }

    let mut version = "3.0".to_string();
    let mut creation_date = String::new();
    let mut originator = String::new();
    let mut segments: Vec<OemSegment> = Vec::new();
    let mut seg: Option<SegBuilder> = None;
    let mut cur_tag: Option<String> = None;
    let mut buf = Vec::new();

    loop {
        match rdr
            .read_event_into(&mut buf)
            .map_err(|e| FormatError::Format(format!("read_oem_xml: xml: {e}")))?
        {
            Event::Eof => break,
            Event::Start(e) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                let local = name.rsplit(':').next().unwrap_or(&name).to_string();
                match local.as_str() {
                    "oem" => {
                        for attr in e.attributes().flatten() {
                            if attr.key.as_ref() == b"version" {
                                version = String::from_utf8_lossy(&attr.value).into_owned();
                            }
                        }
                    }
                    "segment" => seg = Some(SegBuilder::default()),
                    "stateVector" => {
                        if let Some(s) = seg.as_mut() {
                            s.cur = Some(StateBuilder::default());
                        }
                    }
                    _ => {}
                }
                cur_tag = Some(local);
            }
            Event::End(e) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                let local = name.rsplit(':').next().unwrap_or(&name).to_string();
                match local.as_str() {
                    "stateVector" => {
                        if let Some(s) = seg.as_mut() {
                            if let Some(b) = s.cur.take() {
                                let need = |o: Option<f64>, n: &str| -> Result<f64, FormatError> {
                                    o.ok_or_else(|| {
                                        FormatError::Format(format!(
                                            "read_oem_xml: stateVector missing {n}"
                                        ))
                                    })
                                };
                                let epoch_jd = need(b.epoch, "EPOCH")?;
                                let st = OemState::new(
                                    epoch_jd,
                                    [need(b.x, "X")?, need(b.y, "Y")?, need(b.z, "Z")?],
                                    [
                                        need(b.vx, "X_DOT")?,
                                        need(b.vy, "Y_DOT")?,
                                        need(b.vz, "Z_DOT")?,
                                    ],
                                );
                                s.states.push(st);
                            }
                        }
                    }
                    "segment" => {
                        if let Some(s) = seg.take() {
                            let need =
                                |o: Option<String>, n: &str| -> Result<String, FormatError> {
                                    o.ok_or_else(|| {
                                        FormatError::Format(format!(
                                            "read_oem_xml: metadata missing {n}"
                                        ))
                                    })
                                };
                            segments.push(OemSegment {
                                metadata: OemMetadata {
                                    object_id: need(s.object_id, "OBJECT_ID")?,
                                    object_name: need(s.object_name, "OBJECT_NAME")?,
                                    ref_frame: need(s.ref_frame, "REF_FRAME")?,
                                    time_system: need(s.time_system, "TIME_SYSTEM")?,
                                    center_name: need(s.center_name, "CENTER_NAME")?,
                                },
                                states: s.states,
                            });
                        }
                    }
                    _ => {}
                }
                cur_tag = None;
            }
            Event::Text(t) => {
                let txt = t
                    .unescape()
                    .map_err(|e| FormatError::Format(format!("read_oem_xml: xml text: {e}")))?
                    .to_string();
                let tag = match &cur_tag {
                    Some(t) => t.as_str(),
                    None => continue,
                };
                if let Some(s) = seg.as_mut() {
                    if let Some(b) = s.cur.as_mut() {
                        match tag {
                            "EPOCH" => b.epoch = iso8601_to_jd(&txt),
                            "X" => b.x = txt.parse().ok(),
                            "Y" => b.y = txt.parse().ok(),
                            "Z" => b.z = txt.parse().ok(),
                            "X_DOT" => b.vx = txt.parse().ok(),
                            "Y_DOT" => b.vy = txt.parse().ok(),
                            "Z_DOT" => b.vz = txt.parse().ok(),
                            _ => {}
                        }
                    } else {
                        match tag {
                            "OBJECT_ID" => s.object_id = Some(txt),
                            "OBJECT_NAME" => s.object_name = Some(txt),
                            "REF_FRAME" => s.ref_frame = Some(txt),
                            "TIME_SYSTEM" => s.time_system = Some(txt),
                            "CENTER_NAME" => s.center_name = Some(txt),
                            _ => {}
                        }
                    }
                } else {
                    match tag {
                        "CREATION_DATE" => creation_date = txt,
                        "ORIGINATOR" => originator = txt,
                        _ => {}
                    }
                }
            }
            _ => {}
        }
        buf.clear();
    }

    if segments.is_empty() {
        return Err(FormatError::Format(
            "read_oem_xml: file contains no segments".into(),
        ));
    }
    Ok(OemFile {
        version,
        creation_date,
        originator,
        segments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_header_and_data_lines() {
        let states = vec![
            OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0]),
            OemState::new(
                2_451_545.0 + 30.0 / 86_400.0,
                [6999.0, 225.0, 0.0],
                [-0.24, 7.49, 0.0],
            ),
        ];
        let mut buf = Vec::new();
        write_oem(
            &mut buf,
            &OemMetadata {
                object_id: "1900-001A".into(),
                object_name: "TEST-SAT".into(),
                ref_frame: "EME2000".into(),
                time_system: "TT".into(),
                center_name: "EARTH".into(),
            },
            &states,
        )
        .unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("CCSDS_OEM_VERS"));
        assert!(text.contains("META_START"));
        assert!(text.contains("META_STOP"));
        assert!(text.contains("EME2000"));
        assert_eq!(text.lines().filter(|l| l.starts_with("2000-")).count(), 2);
    }

    #[test]
    fn jd_2451545_is_j2000() {
        // J2000 = 2000-01-01 12:00:00 (TT).
        let s = jd_to_iso8601(2_451_545.0);
        assert!(s.starts_with("2000-01-01T12:00"), "{}", s);
    }

    #[test]
    fn iso_jd_round_trip_at_j2000() {
        let jd = 2_451_545.0;
        let iso = jd_to_iso8601(jd);
        let back = iso8601_to_jd(&iso).expect("parse iso");
        assert!((back - jd).abs() < 1e-9, "{back} vs {jd}");
    }

    #[test]
    fn write_then_read_round_trip() {
        let states = vec![
            OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0]),
            OemState::new(
                2_451_545.0 + 30.0 / 86_400.0,
                [6999.0, 225.0, 0.0],
                [-0.24, 7.49, 0.0],
            ),
        ];
        let meta = OemMetadata {
            object_id: "1900-001A".into(),
            object_name: "TEST-SAT".into(),
            ref_frame: "EME2000".into(),
            time_system: "TT".into(),
            center_name: "EARTH".into(),
        };
        let mut buf = Vec::new();
        write_oem(&mut buf, &meta, &states).unwrap();
        let parsed = read_oem(&buf[..]).expect("parse OEM");
        assert_eq!(parsed.version, "3.0");
        assert_eq!(parsed.segments.len(), 1);
        let seg = &parsed.segments[0];
        assert_eq!(seg.metadata.object_id, "1900-001A");
        assert_eq!(seg.metadata.ref_frame, "EME2000");
        assert_eq!(seg.states.len(), 2);
        for (a, b) in seg.states.iter().zip(states.iter()) {
            assert!((a.position_km[0] - b.position_km[0]).abs() < 1e-3);
            assert!((a.position_km[1] - b.position_km[1]).abs() < 1e-3);
            assert!((a.position_km[2] - b.position_km[2]).abs() < 1e-3);
            assert!((a.velocity_km_s[0] - b.velocity_km_s[0]).abs() < 1e-6);
            assert!((a.velocity_km_s[1] - b.velocity_km_s[1]).abs() < 1e-6);
            assert!((a.velocity_km_s[2] - b.velocity_km_s[2]).abs() < 1e-6);
        }
    }

    #[test]
    fn read_oem_rejects_missing_meta_field() {
        let txt = "CCSDS_OEM_VERS = 3.0\n\
                   META_START\n\
                   OBJECT_NAME = X\n\
                   META_STOP\n";
        let err = read_oem(txt.as_bytes()).unwrap_err();
        match err {
            FormatError::Format(msg) => assert!(msg.contains("OBJECT_ID"), "{msg}"),
            other => panic!("expected Format error, got {other:?}"),
        }
    }

    #[test]
    fn xml_round_trip() {
        let states = vec![
            OemState::new(2_451_545.0, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0]),
            OemState::new(
                2_451_545.0 + 30.0 / 86_400.0,
                [6999.0, 225.0, 0.0],
                [-0.24, 7.49, 0.0],
            ),
        ];
        let meta = OemMetadata {
            object_id: "1900-001A".into(),
            object_name: "TEST-SAT".into(),
            ref_frame: "EME2000".into(),
            time_system: "TT".into(),
            center_name: "EARTH".into(),
        };
        let mut buf = Vec::new();
        write_oem_xml(&mut buf, &meta, &states).unwrap();
        let parsed = read_oem_xml(&buf[..]).expect("parse OEM XML");
        assert_eq!(parsed.segments.len(), 1);
        let seg = &parsed.segments[0];
        assert_eq!(seg.metadata.ref_frame, "EME2000");
        assert_eq!(seg.states.len(), 2);
        for (a, b) in seg.states.iter().zip(states.iter()) {
            assert!((a.position_km[0] - b.position_km[0]).abs() < 1e-3);
            assert!((a.velocity_km_s[1] - b.velocity_km_s[1]).abs() < 1e-6);
        }
    }
}
