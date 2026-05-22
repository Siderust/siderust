//! # RINEX observation reader
//!
//! ## Scientific scope
//!
//! RINEX observation files provide the time-tagged code and carrier
//! measurements that drive GNSS POD. This module implements the compact
//! subset needed by the current workspace: a small selection of RINEX 3
//! header fields and dual-frequency observation values for the supported
//! systems.
//!
//! It is designed for deterministic ingestion of the dual-frequency
//! observables needed by the POD pipeline rather than exhaustive
//! archive support. Unsupported observables and header records are
//! skipped instead of being approximated.
//!
//! ## Technical scope
//!
//! The public items are `ObsEpoch`, `RinexObs`, and `read_rinex_obs`.
//! Parsed epochs retain the file-native observation codes and scalar values
//! while associating them with typed epoch information suitable for later
//! observation modelling.
//!
//! Bias handling, ambiguity resolution, and troposphere/ionosphere
//! modelling are beyond the scope of this parser.
//!
//! ## References
//!
//! - International GNSS Service / RTCM. (2020). RINEX: The Receiver
//!   Independent Exchange Format, Version 3.05.
//! - Misra, P., & Enge, P. (2012). Global Positioning System: Signals,
//!   Measurements, and Performance (2nd ed.). Ganga-Jamuna Press.
use super::FormatError;
use chrono::{DateTime, Datelike, NaiveDate, Timelike, Utc as ChronoUtc};
use qtty::length::Meters;
use qtty::time::Seconds;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};
use tempoch::{Time, UTC};

/// One observation epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct ObsEpoch {
    /// UTC epoch.
    pub time: Time<UTC>,
    /// Per-satellite observations: `sat -> {obs_code -> value}`.
    pub satellites: HashMap<String, HashMap<String, f64>>,
}

/// Parsed RINEX OBS file.
#[derive(Debug, Clone)]
pub struct RinexObs {
    /// Marker name (station name).
    pub marker: String,
    /// Approximate XYZ in metres (ITRF).
    pub approx_xyz_m: Option<[Meters; 3]>,
    /// Sampling interval, if declared.
    pub interval_s: Option<Seconds>,
    /// Per-system observation type list.
    pub obs_types: HashMap<char, Vec<String>>,
    /// Epoch records.
    pub epochs: Vec<ObsEpoch>,
}

/// Parse a RINEX OBS file.
pub fn read_rinex_obs<R: Read>(rdr: R) -> Result<RinexObs, FormatError> {
    let mut br = BufReader::new(rdr);
    let mut line = String::new();

    let mut marker = String::new();
    let mut approx_xyz_m: Option<[Meters; 3]> = None;
    let mut interval_s: Option<Seconds> = None;
    let mut obs_types: HashMap<char, Vec<String>> = HashMap::new();

    // Header.
    loop {
        line.clear();
        if br.read_line(&mut line)? == 0 {
            return Err(FormatError::Format("rinex_obs: unexpected EOF".into()));
        }
        let label = line.get(60..).unwrap_or("").trim();
        let body = line.get(..60).unwrap_or("");
        match label {
            "MARKER NAME" => marker = body.trim().to_string(),
            "APPROX POSITION XYZ" => {
                let parts: Vec<f64> = body
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() == 3 {
                    approx_xyz_m = Some([
                        Meters::new(parts[0]),
                        Meters::new(parts[1]),
                        Meters::new(parts[2]),
                    ]);
                }
            }
            "INTERVAL" => {
                interval_s = body
                    .split_whitespace()
                    .next()
                    .and_then(|s| s.parse().ok())
                    .map(Seconds::new);
            }
            "SYS / # / OBS TYPES" => {
                let bytes = body.as_bytes();
                if bytes.is_empty() {
                    continue;
                }
                let sys = bytes[0] as char;
                if !sys.is_ascii_alphabetic() {
                    continue;
                }
                let count: usize = body[1..6].trim().parse().unwrap_or(0);
                let mut types: Vec<String> = body[7..]
                    .split_whitespace()
                    .map(|s| s.to_string())
                    .collect();
                while types.len() < count {
                    let mut cont = String::new();
                    if br.read_line(&mut cont)? == 0 {
                        break;
                    }
                    let c_body = cont.get(..60).unwrap_or("");
                    types.extend(c_body.split_whitespace().map(|s| s.to_string()));
                }
                types.truncate(count);
                obs_types.insert(sys, types);
            }
            "END OF HEADER" => break,
            _ => {}
        }
    }

    // Epochs.
    let mut epochs = Vec::new();
    loop {
        line.clear();
        if br.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end();
        if !trimmed.starts_with('>') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 9 {
            continue;
        }
        let year: i32 = parts[1].parse().unwrap_or(0);
        let month: u32 = parts[2].parse().unwrap_or(0);
        let day: u32 = parts[3].parse().unwrap_or(0);
        let hour: u32 = parts[4].parse().unwrap_or(0);
        let minute: u32 = parts[5].parse().unwrap_or(0);
        let second_f: f64 = parts[6].parse().unwrap_or(0.0);
        let second_int = second_f as u32;
        let nanos = ((second_f - second_int as f64) * 1e9).round() as u32;
        let naive = NaiveDate::from_ymd_opt(year, month, day)
            .and_then(|d| d.and_hms_nano_opt(hour, minute, second_int, nanos))
            .ok_or_else(|| FormatError::Format("rinex_obs: invalid epoch date".into()))?;
        let epoch_time =
            Time::<UTC>::try_from_chrono(DateTime::from_naive_utc_and_offset(naive, ChronoUtc))
                .map_err(|e| {
                    FormatError::Format(format!("rinex_obs: epoch UTC conversion: {e}"))
                })?;
        let n_sat: usize = parts[8].parse().unwrap_or(0);

        let mut sats = HashMap::new();
        for _ in 0..n_sat {
            line.clear();
            if br.read_line(&mut line)? == 0 {
                break;
            }
            let row = line.trim_end();
            if row.len() < 3 {
                continue;
            }
            let sat = row[..3].to_string();
            let sys = sat.chars().next().unwrap_or(' ');
            let types = match obs_types.get(&sys) {
                Some(t) => t.clone(),
                None => Vec::new(),
            };
            let mut vals: HashMap<String, f64> = HashMap::new();
            // RINEX 3: each observable occupies 16 columns (F14.3 + 2 status
            // chars). We slice that and parse loosely.
            let payload = &row[3..];
            for (i, ty) in types.iter().enumerate() {
                let start = i * 16;
                if start >= payload.len() {
                    break;
                }
                let end = (start + 14).min(payload.len());
                let raw = payload[start..end].trim();
                if raw.is_empty() {
                    continue;
                }
                if let Ok(v) = raw.parse::<f64>() {
                    vals.insert(ty.clone(), v);
                }
            }
            sats.insert(sat, vals);
        }
        epochs.push(ObsEpoch {
            time: epoch_time,
            satellites: sats,
        });
    }

    Ok(RinexObs {
        marker,
        approx_xyz_m,
        interval_s,
        obs_types,
        epochs,
    })
}

/// Write a RINEX 3.04 observation file matching the layout produced by
/// [`read_rinex_obs`].
///
/// The writer emits a minimal but spec-conformant RINEX 3 header (VERSION /
/// TYPE, MARKER NAME, optional APPROX POSITION XYZ, SYS / # / OBS TYPES per
/// system, optional INTERVAL, END OF HEADER) followed by one epoch block
/// per [`ObsEpoch`]. Per-satellite observations are written in the order
/// declared by `obs_types[&sys]`; missing observations are written as
/// blanks (the standard sentinel for "not observed").
///
/// Round-trip parity is **canonicalised**: header line order is normalised,
/// trailing whitespace is stripped, and observation columns are emitted in
/// the standard 16-column F14.3 layout. Numerical values round-trip to
/// 1e-3 m (code) / 1e-3 cycles (carrier), which exceeds RINEX 3 native
/// precision.
///
/// # Errors
///
/// Returns [`FormatError::Io`] on write failures and [`FormatError::Format`]
/// if an epoch cannot be projected onto the proleptic Gregorian calendar.
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::obs::{read_rinex_obs, write_rinex_obs};
///
/// let src = concat!(
///     "     3.04           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE\n",
///     "TINY                                                        MARKER NAME\n",
///     "G    1 C1C                                                  SYS / # / OBS TYPES\n",
///     "                                                            END OF HEADER\n",
///     "> 2024 01 01 00 00  0.0000000  0  1\n",
///     "G01  20100000.000\n",
/// );
/// let obs = read_rinex_obs(src.as_bytes()).unwrap();
/// let mut buf = Vec::new();
/// write_rinex_obs(&mut buf, &obs).unwrap();
/// let parsed = read_rinex_obs(&buf[..]).unwrap();
/// assert_eq!(parsed.marker, obs.marker);
/// assert_eq!(parsed.epochs.len(), obs.epochs.len());
/// ```
pub fn write_rinex_obs<W: Write>(w: &mut W, obs: &RinexObs) -> Result<(), FormatError> {
    fn header_line(w: &mut impl Write, body: &str, label: &str) -> std::io::Result<()> {
        let body = if body.len() > 60 { &body[..60] } else { body };
        writeln!(w, "{:<60}{}", body, label)
    }
    header_line(
        w,
        "     3.04           OBSERVATION DATA    M (MIXED)",
        "RINEX VERSION / TYPE",
    )?;
    header_line(w, &obs.marker, "MARKER NAME")?;
    if let Some(xyz) = &obs.approx_xyz_m {
        let body = format!(
            "{:14.4}{:14.4}{:14.4}",
            xyz[0].value(),
            xyz[1].value(),
            xyz[2].value()
        );
        header_line(w, &body, "APPROX POSITION XYZ")?;
    }
    let mut systems: Vec<&char> = obs.obs_types.keys().collect();
    systems.sort();
    for sys in systems {
        let types = &obs.obs_types[sys];
        let mut body = format!("{}{:>5}", sys, types.len());
        for (i, t) in types.iter().enumerate() {
            // RINEX 3: 13 obs codes per header line, 4 chars each.
            if i > 0 && i % 13 == 0 {
                header_line(w, &body, "SYS / # / OBS TYPES")?;
                body = format!("{:>6}", "");
            }
            body.push_str(&format!(" {:<3}", t));
        }
        header_line(w, &body, "SYS / # / OBS TYPES")?;
    }
    if let Some(interval) = obs.interval_s {
        let body = format!("{:10.3}", interval.value());
        header_line(w, &body, "INTERVAL")?;
    }
    header_line(w, "", "END OF HEADER")?;

    for ep in &obs.epochs {
        let dt = ep
            .time
            .try_to_chrono()
            .map_err(|e| FormatError::Format(format!("rinex_obs: epoch to chrono failed: {e}")))?;
        let second = dt.second() as f64 + dt.nanosecond() as f64 / 1e9;
        writeln!(
            w,
            "> {:4} {:02} {:02} {:02} {:02} {:10.7}  0 {:>2}",
            dt.year(),
            dt.month(),
            dt.day(),
            dt.hour(),
            dt.minute(),
            second,
            ep.satellites.len(),
        )?;
        let mut sats: Vec<&String> = ep.satellites.keys().collect();
        sats.sort();
        for sat in sats {
            let vals = &ep.satellites[sat];
            let sys = sat.chars().next().unwrap_or(' ');
            let types = obs.obs_types.get(&sys).cloned().unwrap_or_default();
            let mut row = format!("{:<3}", sat);
            for ty in &types {
                match vals.get(ty) {
                    Some(v) => row.push_str(&format!("{:14.3}  ", v)),
                    None => row.push_str(&" ".repeat(16)),
                }
            }
            writeln!(w, "{}", row.trim_end())?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
     3.04           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE
TEST-MARK                                                   MARKER NAME
  1234567.8901  -2345678.9012   5678901.2345                APPROX POSITION XYZ
G    2 C1C L1C                                              SYS / # / OBS TYPES
    30                                                      INTERVAL
                                                            END OF HEADER
> 2024 01 01 00 00  0.0000000  0  1
G01  20123456.789         123456789.012
> 2024 01 01 00 00 30.0000000  0  1
G01  20124000.000         123457000.000
";

    #[test]
    fn parses_header_and_epochs() {
        let r = read_rinex_obs(SAMPLE.as_bytes()).unwrap();
        assert_eq!(r.marker, "TEST-MARK");
        assert_eq!(r.interval_s.map(|s| s.value()), Some(30.0));
        assert_eq!(r.obs_types[&'G'], vec!["C1C", "L1C"]);
        assert_eq!(r.epochs.len(), 2);
        let e0 = &r.epochs[0];
        assert!(e0.satellites.contains_key("G01"));
        let g01 = &e0.satellites["G01"];
        assert!((g01["C1C"] - 20_123_456.789).abs() < 1e-3);
    }

    #[test]
    fn round_trips_through_writer() {
        let r = read_rinex_obs(SAMPLE.as_bytes()).unwrap();
        let mut buf = Vec::new();
        write_rinex_obs(&mut buf, &r).unwrap();
        let r2 = read_rinex_obs(&buf[..]).unwrap();
        assert_eq!(r2.marker, r.marker);
        assert_eq!(r2.epochs.len(), r.epochs.len());
        assert_eq!(r2.obs_types, r.obs_types);
        assert_eq!(
            r2.interval_s.map(|s| s.value()),
            r.interval_s.map(|s| s.value())
        );
        for (a, b) in r.epochs.iter().zip(r2.epochs.iter()) {
            assert_eq!(
                a.satellites
                    .keys()
                    .collect::<std::collections::BTreeSet<_>>(),
                b.satellites
                    .keys()
                    .collect::<std::collections::BTreeSet<_>>()
            );
            for (sat, vals_a) in &a.satellites {
                let vals_b = &b.satellites[sat];
                for (k, v) in vals_a {
                    assert!(
                        (vals_b[k] - v).abs() < 1e-3,
                        "{sat}/{k} {v} vs {}",
                        vals_b[k]
                    );
                }
            }
        }
    }
}
