//! # CRD SLR observation reader
//!
//! ## Scientific scope
//!
//! The Consolidated Laser Ranging Data (CRD) format is the primary interchange
//! format for Satellite Laser Ranging (SLR) observations published by the
//! International Laser Ranging Service (ILRS). This module supports CRD
//! versions 1 and 2 with both strict and permissive parse modes.
//!
//! The primary POD product is the **normal-point** record (type `11`), which
//! is a compressed, time-averaged two-way range observation. Full-rate records
//! (type `10`) are also parsed. Station and target metadata (H2/H3 header
//! blocks) are captured in typed structs.
//!
//! Station motion, atmospheric delay, range residual analysis, and any
//! corrections beyond what is in the file are outside this module.
//!
//! ## Technical scope
//!
//! Public entry points are [`read_crd`], [`parse_crd`], and
//! [`parse_crd_with_mode`]. All three return a [`CrdFile`] containing:
//!
//! - [`CrdStation`] / [`CrdTarget`] metadata from H2/H3.
//! - A [`Vec<NormalPoint>`] of typed normal-point records.
//! - A [`Vec<CrdRange>`] of raw range records (both types 10 and 11) for
//!   backward compatibility.
//!
//! The one-way slant range in [`NormalPoint::range_m`] is derived as
//! `c × TOF / 2` where `c = 299 792 458 m s⁻¹`.
//!
//! ## References
//!
//! - International Laser Ranging Service. (2022). Consolidated Laser
//!   Ranging Data Format Specification (v2.01).
//!   <https://ilrs.gsfc.nasa.gov/docs/2022/ILRS_CRD_Format_v2.01.pdf>
//! - Pearlman, M. R., Noll, C. E., et al. (2019). The ILRS: Current status
//!   and future prospects. Journal of Geodesy, 93, 2161–2180.
use super::{FileLocation, FormatError, ParseMode};
use chrono::{DateTime, NaiveDate, Utc as ChronoUtc};
use qtty::length::Meters;
use qtty::time::Seconds;
use std::fs;
use std::path::Path;
use tempoch::{Time, UTC};

/// Speed of light in m/s (IAU 2012 definition).
const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

// ── Station / Target metadata ─────────────────────────────────────────────────

/// Station metadata from the CRD H2 header record.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::{CrdStation, parse_crd};
/// let txt = "H2 GRAZ 7839 1 1 0\nH8\n";
/// let f = parse_crd(txt).unwrap();
/// assert_eq!(f.station.name, "GRAZ");
/// assert_eq!(f.station.cdp_pad, 7839);
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CrdStation {
    /// Station name (e.g. `"GRAZ"`).
    pub name: String,
    /// CDP Pad number.
    pub cdp_pad: i32,
    /// System number.
    pub sys_no: i32,
    /// Occupancy sequence number.
    pub occ_no: i32,
    /// Station time zone (hours from UTC, usually 0).
    pub time_zone: i32,
}

/// Target (satellite) metadata from the CRD H3 header record.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::{CrdTarget, parse_crd};
/// let txt = "H3 lageos1 1155 7603901 0 1\nH8\n";
/// let f = parse_crd(txt).unwrap();
/// assert_eq!(f.target.name, "lageos1");
/// assert_eq!(f.target.sic, 1155);
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CrdTarget {
    /// Target name (e.g. `"lageos1"`).
    pub name: String,
    /// Satellite Identification Code.
    pub sic: i32,
    /// COSPAR/NORAD-style identifier.
    pub norad: String,
    /// Spacecraft flag (0 = passive, 1 = active).
    pub sc_flag: i32,
    /// Epoch identifier.
    pub epoch_id: i32,
}

// ── Typed normal-point record ─────────────────────────────────────────────────

/// A typed normal-point observation (CRD record type `11`).
///
/// A normal point compresses many raw returns within a time bin into a single
/// smoothed two-way range measurement. This is the primary POD input product
/// for SLR stations.
///
/// The one-way slant range [`range_m`](NormalPoint::range_m) is computed as
/// `c × time_of_flight / 2`.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::parse_crd;
///
/// let txt = "\
/// H1 CRD 2 2024 01 01 00\n\
/// H2 GRAZ 7839 1 1 0\n\
/// H3 lageos1 1155 7603901 0 1\n\
/// H4 0 2024 01 01 08 00 00 12 00 00 0 0 0 0 1 0\n\
/// C0 0 std 1 5 0 0 1 0 0 0 0 0 0 0\n\
/// 11 28800.0 0.05123456789 std 0 0 1 0 5 8.0 0.1 0.0 0.0 90.0\n\
/// H8\n";
/// let f = parse_crd(txt).unwrap();
/// assert_eq!(f.normal_points.len(), 1);
/// let np = &f.normal_points[0];
/// assert!((np.time_of_flight.value() - 0.05123456789).abs() < 1e-12);
/// // One-way range ≈ 7 678 km
/// assert!(np.range_m.value() > 7_000_000.0);
/// ```
#[derive(Debug, Clone)]
pub struct NormalPoint {
    /// Seconds of day (UTC) at the observation epoch.
    pub seconds_of_day: Seconds,
    /// UTC epoch (session date midnight + SOD). `None` if H4 was absent.
    pub epoch: Option<Time<UTC>>,
    /// Two-way time-of-flight in seconds.
    pub time_of_flight: Seconds,
    /// One-way slant range in metres (`= c × TOF / 2`).
    pub range_m: Meters,
    /// System configuration ID from the preceding C0 record.
    pub system_config_id: String,
    /// Epoch event (0 = centre, 1 = start, 2 = end of bin).
    pub epoch_event: u8,
    /// Filter / calibration flag.
    pub filter_flag: u8,
    /// Data quality indicator.
    pub data_quality: u8,
    /// Format flag.
    pub format_flag: u8,
    /// Number of raw ranges contributing to this normal point.
    pub num_raws: Option<u32>,
    /// RMS of raw two-way ranges within the bin (m).
    pub bin_rms_m: Option<Meters>,
    /// Normal-point bin size (s).
    pub bin_size_s: Option<Seconds>,
    /// Return rate within the bin (%).
    pub return_rate: Option<f64>,
}

// ── Backward-compatible raw range record ──────────────────────────────────────

/// A raw SLR range record from CRD type `10` (full-rate) or `11` (normal-point).
///
/// Kept for backward compatibility. New code should prefer [`NormalPoint`]
/// for type-11 records.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::parse_crd;
/// let txt = "H2 S 0 0 0 0\nH3 t 0 0 0 1\nH4 0 2024 1 1 0 0 0 0 0 0 0 0 0 0 1 0\n\
///            11 300.0 0.08 std 0 0 1 0\nH8\n";
/// let f = parse_crd(txt).unwrap();
/// assert_eq!(f.ranges[0].record_type, 11);
/// ```
#[derive(Debug, Clone)]
pub struct CrdRange {
    /// Seconds of day (UTC) at the epoch of the range observation.
    pub seconds_of_day: Seconds,
    /// Two-way time-of-flight.
    pub time_of_flight: Seconds,
    /// System configuration ID.
    pub system_config_id: String,
    /// `10` for full-rate, `11` for normal-point.
    pub record_type: u8,
}

// ── CrdFile ───────────────────────────────────────────────────────────────────

/// Parsed CRD file (header + observations).
///
/// Contains both the high-level typed view ([`normal_points`](CrdFile::normal_points))
/// and the backward-compatible raw [`ranges`](CrdFile::ranges).
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::parse_crd;
/// let txt = "H2 GRAZ 7839 1 1 0\nH3 lageos1 1155 7603901 0 1\nH8\n";
/// let f = parse_crd(txt).unwrap();
/// assert_eq!(f.station_name, "GRAZ");
/// assert_eq!(f.target.sic, 1155);
/// ```
#[derive(Debug, Clone, Default)]
pub struct CrdFile {
    /// Station name — mirrors `station.name` for backward compatibility.
    pub station_name: String,
    /// CDP/Pad number — mirrors `station.cdp_pad`.
    pub station_cdp_pad: i32,
    /// Satellite name — mirrors `target.name`.
    pub satellite_name: String,
    /// SIC — mirrors `target.sic`.
    pub satellite_sic: i32,
    /// NORAD/COSPAR ID — mirrors `target.norad`.
    pub satellite_norad: String,
    /// Session start (UTC midnight) from the H4 record.
    pub session_date: Option<Time<UTC>>,
    /// Raw range records (types 10 and 11) — backward-compatible view.
    pub ranges: Vec<CrdRange>,
    /// CRD format version string from H1 (e.g. `"2"`).
    pub format_version: String,
    /// Typed station metadata.
    pub station: CrdStation,
    /// Typed target (satellite) metadata.
    pub target: CrdTarget,
    /// Normal-point records (type 11) — primary POD product.
    pub normal_points: Vec<NormalPoint>,
    /// Parse mode used.
    pub parse_mode: ParseMode,
}

// ── Public entry points ───────────────────────────────────────────────────────

/// Read a CRD file from disk using [`ParseMode::Strict`].
///
/// # Examples
///
/// ```no_run
/// use siderust::formats::ilrs::crd::read_crd;
/// let f = read_crd("station.crd").unwrap();
/// println!("{} normal points", f.normal_points.len());
/// ```
pub fn read_crd<P: AsRef<Path>>(path: P) -> Result<CrdFile, FormatError> {
    let text = fs::read_to_string(path.as_ref())?;
    parse_crd_impl(&text, ParseMode::Strict, Some(path.as_ref().into()))
}

/// Parse a CRD string using [`ParseMode::Strict`].
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::crd::parse_crd;
/// let f = parse_crd("H2 GRAZ 7839 1 1 0\nH8\n").unwrap();
/// assert_eq!(f.station_name, "GRAZ");
/// ```
pub fn parse_crd(text: &str) -> Result<CrdFile, FormatError> {
    parse_crd_impl(text, ParseMode::Strict, None)
}

/// Parse a CRD string with an explicit [`ParseMode`].
///
/// In `Strict` mode, a malformed range record (missing SOD or TOF) is a hard
/// error. In `Permissive` mode, such records are silently skipped.
///
/// # Examples
///
/// ```
/// use siderust::formats::{ParseMode, ilrs::crd::parse_crd_with_mode};
/// // Permissive: a broken "11" record is skipped rather than failing.
/// let txt = "H2 GRAZ 7839 1 1 0\nH3 lageos1 1155 7603901 0 1\n11 BROKEN\nH8\n";
/// let f = parse_crd_with_mode(txt, ParseMode::Permissive).unwrap();
/// assert_eq!(f.normal_points.len(), 0); // skipped
/// ```
pub fn parse_crd_with_mode(text: &str, mode: ParseMode) -> Result<CrdFile, FormatError> {
    parse_crd_impl(text, mode, None)
}

// ── Internal parser ───────────────────────────────────────────────────────────

fn parse_crd_impl(
    text: &str,
    mode: ParseMode,
    path: Option<std::path::PathBuf>,
) -> Result<CrdFile, FormatError> {
    let mut out = CrdFile {
        parse_mode: mode,
        ..Default::default()
    };
    let mut current_sys = String::from("std");
    let mut h4_date: Option<(i32, u32, u32)> = None;

    for (line_idx, raw) in text.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = raw.trim_end();
        if line.is_empty() {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let tag = match tokens.next() {
            Some(t) => t,
            None => continue,
        };

        match tag {
            "H1" => {
                // H1 CRD <version> <year> <month> <day> <hour>
                let _fmt = tokens.next(); // "CRD"
                if let Some(v) = tokens.next() {
                    out.format_version = v.to_string();
                }
            }
            "H2" => {
                // H2 <station_name> <cdp_pad> <sys_no> <occ_no> <time_zone>
                let name = tokens.next().unwrap_or("").to_string();
                let cdp_pad = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let sys_no = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let occ_no = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let time_zone = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                out.station = CrdStation {
                    name: name.clone(),
                    cdp_pad,
                    sys_no,
                    occ_no,
                    time_zone,
                };
                out.station_name = name;
                out.station_cdp_pad = cdp_pad;
            }
            "H3" => {
                // H3 <target_name> <sic> <norad> <sc_flag> <epoch_id>
                let name = tokens.next().unwrap_or("").to_string();
                let sic: i32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let norad = tokens.next().unwrap_or("").to_string();
                let sc_flag: i32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let epoch_id: i32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                out.target = CrdTarget {
                    name: name.clone(),
                    sic,
                    norad: norad.clone(),
                    sc_flag,
                    epoch_id,
                };
                out.satellite_name = name;
                out.satellite_sic = sic;
                out.satellite_norad = norad;
            }
            "H4" => {
                // H4 <type> <year> <month> <day> <hour> <min> <sec> ...
                let _kind = tokens.next();
                let year: i32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let month: u32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let day: u32 = tokens.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                if let Some(naive) =
                    NaiveDate::from_ymd_opt(year, month, day).and_then(|d| d.and_hms_opt(0, 0, 0))
                {
                    let dt = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
                    if let Ok(t) = Time::<UTC>::try_from_chrono(dt) {
                        out.session_date = Some(t);
                        h4_date = Some((year, month, day));
                    }
                }
            }
            "H5" | "H8" | "H9" => { /* session/file end markers — ignored */ }
            "C0" => {
                // C0 <detail_type> <config_id> ...
                let _detail = tokens.next();
                if let Some(s) = tokens.next() {
                    current_sys = s.to_string();
                }
            }
            "C1" | "C2" | "C3" | "C4" => { /* other config records — skipped */ }
            "10" | "11" => {
                let record_type: u8 = if tag == "11" { 11 } else { 10 };

                let sod_str = tokens.next();
                let tof_str = tokens.next();

                let (sod, tof) = match (
                    sod_str.and_then(|s| s.parse::<f64>().ok()),
                    tof_str.and_then(|s| s.parse::<f64>().ok()),
                ) {
                    (Some(s), Some(t)) => (s, t),
                    (None, _) => {
                        let loc = FileLocation::new(path.clone(), Some(line_no), None);
                        let err = FormatError::located(
                            "CRD v2 §4.1",
                            loc,
                            format!("record {tag}: missing seconds-of-day"),
                        );
                        if mode == ParseMode::Strict {
                            return Err(err);
                        }
                        continue;
                    }
                    (_, None) => {
                        let loc = FileLocation::new(path.clone(), Some(line_no), None);
                        let err = FormatError::located(
                            "CRD v2 §4.1",
                            loc,
                            format!("record {tag}: missing time-of-flight"),
                        );
                        if mode == ParseMode::Strict {
                            return Err(err);
                        }
                        continue;
                    }
                };

                // Gather remaining optional fields for normal points.
                let rest: Vec<&str> = tokens.collect();

                // Backward-compatible raw range record.
                out.ranges.push(CrdRange {
                    seconds_of_day: Seconds::new(sod),
                    time_of_flight: Seconds::new(tof),
                    system_config_id: current_sys.clone(),
                    record_type,
                });

                if record_type == 11 {
                    // Build epoch directly from H4 date + SOD components to
                    // avoid floating-point precision loss in the JD round-trip.
                    let epoch = h4_date.and_then(|(y, mo, d)| {
                        let whole = sod as i64;
                        let subsec_nanos = ((sod - whole as f64) * 1_000_000_000.0).round() as u32;
                        let hh = (whole / 3600) as u32;
                        let mm = ((whole % 3600) / 60) as u32;
                        let ss = (whole % 60) as u32;
                        let naive = NaiveDate::from_ymd_opt(y, mo, d)
                            .and_then(|nd| nd.and_hms_nano_opt(hh, mm, ss, subsec_nanos))?;
                        let dt = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
                        Time::<UTC>::try_from_chrono(dt).ok()
                    });
                    let range_m = Meters::new(tof * SPEED_OF_LIGHT_M_S / 2.0);

                    // Parse optional fields from the "11" record.
                    // CRD v2 §4.1 record-11 field order (after TOF):
                    // sys_cfg_id epoch_event filter_flag data_quality
                    // format_flag num_raws bin_rms(ps) bin_skew bin_kurtosis
                    // bin_peak return_rate [detector_channel]
                    let sys_id = rest
                        .first()
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| current_sys.clone());
                    let epoch_event: u8 = rest.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let filter_flag: u8 = rest.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let data_quality: u8 = rest.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let format_flag: u8 = rest.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let num_raws: Option<u32> = rest.get(5).and_then(|s| s.parse().ok());
                    // bin_rms is in picoseconds in the file; convert to metres.
                    let bin_rms_m: Option<Meters> = rest
                        .get(6)
                        .and_then(|s| s.parse::<f64>().ok())
                        .map(|ps| Meters::new(ps * 1e-12 * SPEED_OF_LIGHT_M_S / 2.0));
                    let return_rate: Option<f64> = rest.get(9).and_then(|s| s.parse().ok());

                    out.normal_points.push(NormalPoint {
                        seconds_of_day: Seconds::new(sod),
                        epoch,
                        time_of_flight: Seconds::new(tof),
                        range_m,
                        system_config_id: sys_id,
                        epoch_event,
                        filter_flag,
                        data_quality,
                        format_flag,
                        num_raws,
                        bin_rms_m,
                        bin_size_s: None,
                        return_rate,
                    });
                }
            }
            // Records 12–52: supplements/met/pointing/cal/stats — skip.
            "12" | "20" | "21" | "30" | "40" | "41" | "50" | "60" => {}
            _ => { /* unknown tag — silently skip in both modes */ }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::ParseMode;
    use super::*;

    #[test]
    fn parses_minimal_crd() {
        let txt = "\
H1 CRD 2 2024 01 01 00\n\
H2 7090 7090 1 01 0\n\
H3 lageos1 1155 7603901 8820 1\n\
H4 0 2024 01 01 00 00 00 00 00 00 0 0 0 0 1 0\n\
C0 0 std\n\
11 12345.000 0.05123456789 std 1 1 1 1 0 0\n\
11 12375.000 0.05123446712 std 1 1 1 1 0 0\n\
H8\n";
        let f = parse_crd(txt).expect("parse");
        assert_eq!(f.station_name, "7090");
        assert_eq!(f.satellite_name, "lageos1");
        use chrono::Datelike;
        let d = f.session_date.unwrap().try_to_chrono().unwrap();
        assert_eq!(d.year(), 2024);
        assert_eq!(d.month(), 1);
        assert_eq!(d.day(), 1);
        assert_eq!(f.ranges.len(), 2);
        assert!((f.ranges[0].time_of_flight.value() - 0.05123456789).abs() < 1e-15);
        assert_eq!(f.ranges[0].record_type, 11);
    }

    #[test]
    fn normal_points_populated() {
        let txt = "\
H1 CRD 2 2024 01 01 00\n\
H2 GRAZ 7839 1 1 0\n\
H3 lageos1 1155 7603901 0 1\n\
H4 0 2024 01 01 08 00 00 12 00 00 0 0 0 0 1 0\n\
C0 0 std 1 5 0 0 1 0 0 0 0 0 0 0\n\
11 28800.0 0.051234567890 std 0 0 1 0 5 8.0 0.0 0.0 0.0 90.0\n\
11 28830.0 0.051198765432 std 0 0 1 0 5 7.5 0.0 0.0 0.0 88.0\n\
H8\n";
        let f = parse_crd(txt).expect("parse");
        assert_eq!(f.normal_points.len(), 2);
        let np = &f.normal_points[0];
        assert!((np.time_of_flight.value() - 0.051234567890).abs() < 1e-12);
        // Range ≈ c * TOF / 2
        let expected_m = 0.051234567890 * SPEED_OF_LIGHT_M_S / 2.0;
        assert!((np.range_m.value() - expected_m).abs() < 1.0);
        assert_eq!(np.num_raws, Some(5));
    }

    #[test]
    fn epoch_computed_from_session_date() {
        let txt = "\
H1 CRD 2 2024 01 01 00\n\
H2 GRAZ 7839 1 1 0\n\
H3 lageos1 1155 7603901 0 1\n\
H4 0 2024 01 01 08 00 00 12 00 00 0 0 0 0 1 0\n\
11 3600.0 0.051 std 0 0 0 0\n\
H8\n";
        let f = parse_crd(txt).expect("parse");
        let np = &f.normal_points[0];
        assert!(np.epoch.is_some(), "epoch should be computed from H4 date");
        use chrono::Timelike;
        let dt = np.epoch.unwrap().try_to_chrono().unwrap();
        assert_eq!(dt.hour(), 1); // 3600 s after midnight = 01:00
    }

    #[test]
    fn station_and_target_populated() {
        let txt = "\
H2 GRAZ 7839 1 2 0\n\
H3 lageos1 1155 7603901 0 1\n\
H8\n";
        let f = parse_crd(txt).expect("parse");
        assert_eq!(f.station.name, "GRAZ");
        assert_eq!(f.station.cdp_pad, 7839);
        assert_eq!(f.station.sys_no, 1);
        assert_eq!(f.station.occ_no, 2);
        assert_eq!(f.target.name, "lageos1");
        assert_eq!(f.target.sic, 1155);
        assert_eq!(f.target.norad, "7603901");
    }

    #[test]
    fn format_version_parsed() {
        let txt = "H1 CRD 2 2024 01 01 00\nH8\n";
        let f = parse_crd(txt).expect("parse");
        assert_eq!(f.format_version, "2");
    }

    #[test]
    fn strict_missing_sod_fails() {
        let txt = "H2 S 0 0 0 0\n11\nH8\n";
        let err = parse_crd(txt).expect_err("strict: missing SOD");
        assert!(matches!(err, FormatError::Located { .. }));
    }

    #[test]
    fn permissive_missing_sod_skips() {
        let txt = "H2 S 0 0 0 0\n11\nH8\n";
        let f = parse_crd_with_mode(txt, ParseMode::Permissive).expect("permissive ok");
        assert_eq!(f.normal_points.len(), 0);
    }

    #[test]
    fn strict_missing_tof_fails() {
        let txt = "H2 S 0 0 0 0\n11 123.0\nH8\n";
        let err = parse_crd(txt).expect_err("strict: missing TOF");
        assert!(matches!(err, FormatError::Located { .. }));
    }

    #[test]
    fn permissive_missing_tof_skips() {
        let txt = "H2 S 0 0 0 0\n11 123.0\nH8\n";
        let f = parse_crd_with_mode(txt, ParseMode::Permissive).expect("permissive ok");
        assert_eq!(f.normal_points.len(), 0);
    }

    #[test]
    fn full_rate_records_parsed() {
        let txt = "H2 S 0 0 0 0\n10 100.0 0.08 std\nH8\n";
        let f = parse_crd(txt).expect("parse full-rate");
        assert_eq!(f.ranges.len(), 1);
        assert_eq!(f.ranges[0].record_type, 10);
        assert_eq!(f.normal_points.len(), 0); // type-10 not added to NPs
    }
}
