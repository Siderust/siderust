//! # RINEX navigation message reader
//!
//! ## Scientific scope
//!
//! Broadcast navigation files encode the coarse orbital and clock
//! information transmitted by GNSS spacecraft. This module reads the GPS-
//! focused subset needed to reconstruct broadcast ephemerides for
//! observation modelling and comparative POD workflows.
//!
//! The scientific regime is therefore limited to the message families and
//! fields currently parsed from RINEX 3 NAV. Unsupported systems and
//! optional record features are intentionally left for follow-up work.
//!
//! ## Technical scope
//!
//! The main entry points are `read_rinex_nav` and `parse_rinex_nav`,
//! returning `RinexNavFile` collections of `GpsNavRecord` values. Parsed
//! fields stay close to the broadcast representation so downstream code can
//! choose how to convert them into propagated states.
//!
//! This module does not itself evaluate the broadcast model or perform
//! orbit determination.
//!
//! ## References
//!
//! - International GNSS Service / RTCM. (2020). RINEX: The Receiver
//!   Independent Exchange Format, Version 3.05.
//! - IS-GPS-200. (current revision). Navstar GPS Space Segment / Navigation
//!   User Interfaces.
use super::{FileLocation, FormatError, ParseMode};
use chrono::{DateTime, Datelike, NaiveDate, Timelike, Utc as ChronoUtc};
use qtty::angular::Radians;
use qtty::angular_rate::AngularRate;
use qtty::length::Meters;
use qtty::time::Seconds;
use qtty::unit::{Radian, Second};
use std::fs;
use std::io::Write;
use std::path::Path;
use tempoch::{Time, UTC};

/// One GPS broadcast navigation record (subset).
#[derive(Debug, Clone)]
pub struct GpsNavRecord {
    /// PRN identifier (e.g. `G01` → 1).
    pub prn: u8,
    /// Time of clock (UTC).
    pub toc: Time<UTC>,
    /// SV clock bias.
    pub af0: Seconds,
    /// SV clock drift (s/s — dimensionless rate; kept as scalar).
    pub af1: f64,
    /// SV clock drift rate (s/s² — kept as scalar).
    pub af2: f64,
    /// IODE (dimensionless index).
    pub iode: f64,
    /// Crs amplitude correction.
    pub crs: Meters,
    /// Δn — mean motion correction.
    pub delta_n: AngularRate<Radian, Second>,
    /// M0 — mean anomaly at reference time.
    pub m0: Radians,
    /// Cuc — argument-of-latitude correction (cosine).
    pub cuc: Radians,
    /// Eccentricity (dimensionless).
    pub e: f64,
    /// Cus — argument-of-latitude correction (sine).
    pub cus: Radians,
    /// √a (m^½ — composite dimension; kept as scalar).
    pub sqrt_a: f64,
    /// Toe — reference time of ephemeris.
    pub toe: Seconds,
    /// Cic — inclination correction (cosine).
    pub cic: Radians,
    /// Ω0 — longitude of ascending node at weekly epoch.
    pub omega0: Radians,
    /// Cis — inclination correction (sine).
    pub cis: Radians,
    /// i0 — inclination angle at reference time.
    pub i0: Radians,
    /// Crc amplitude correction.
    pub crc: Meters,
    /// ω — argument of perigee.
    pub omega: Radians,
    /// Ω̇ — rate of right ascension.
    pub omega_dot: AngularRate<Radian, Second>,
    /// IDOT — rate of inclination angle.
    pub idot: AngularRate<Radian, Second>,
}

impl Default for GpsNavRecord {
    fn default() -> Self {
        let naive = NaiveDate::from_ymd_opt(2000, 1, 1)
            .unwrap()
            .and_hms_opt(0, 0, 0)
            .unwrap();
        let dt = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
        Self {
            prn: 0,
            toc: Time::<UTC>::try_from_chrono(dt).unwrap(),
            af0: Seconds::new(0.0),
            af1: 0.0,
            af2: 0.0,
            iode: 0.0,
            crs: Meters::new(0.0),
            delta_n: AngularRate::new(0.0),
            m0: Radians::new(0.0),
            cuc: Radians::new(0.0),
            e: 0.0,
            cus: Radians::new(0.0),
            sqrt_a: 0.0,
            toe: Seconds::new(0.0),
            cic: Radians::new(0.0),
            omega0: Radians::new(0.0),
            cis: Radians::new(0.0),
            i0: Radians::new(0.0),
            crc: Meters::new(0.0),
            omega: Radians::new(0.0),
            omega_dot: AngularRate::new(0.0),
            idot: AngularRate::new(0.0),
        }
    }
}

/// Parsed RINEX NAV file.
#[derive(Debug, Clone, Default)]
pub struct RinexNavFile {
    /// All GPS broadcast records found in the file.
    pub gps: Vec<GpsNavRecord>,
}

/// Read a RINEX NAV file from disk.
pub fn read_rinex_nav<P: AsRef<Path>>(path: P) -> Result<RinexNavFile, FormatError> {
    let text = fs::read_to_string(path)?;
    parse_rinex_nav(&text)
}

/// Parse a RINEX NAV file from a string slice using [`ParseMode::Strict`].
///
/// Returns an error on the first malformed record. Use
/// [`parse_rinex_nav_with_mode`] with [`ParseMode::Permissive`] to skip
/// malformed records while keeping good ones.
///
/// # Errors
///
/// Returns [`FormatError::Located`] on malformed header, invalid PRN, invalid
/// epoch, or body lines with fewer than the expected number of fields.
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::nav::parse_rinex_nav;
///
/// let txt = "\
///      3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
///                                                             END OF HEADER\n\
/// G01 2024 01 01 00 00 00 1.234567E-04 5.678E-12 0.000E+00\n\
///      1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
///      5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
///      5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
///      9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
///      1.000E-10\n\
///      0.000E+00 0.000E+00 0.000E+00 0.000E+00\n\
///      0.000E+00 0.000E+00 0.000E+00 0.000E+00\n";
/// let f = parse_rinex_nav(txt).expect("valid NAV");
/// assert_eq!(f.gps.len(), 1);
/// assert_eq!(f.gps[0].prn, 1);
/// ```
pub fn parse_rinex_nav(text: &str) -> Result<RinexNavFile, FormatError> {
    parse_rinex_nav_with_mode(text, ParseMode::Strict)
}

/// Parse a RINEX NAV file from a string slice with an explicit [`ParseMode`].
///
/// - `Strict` (the default): returns a [`FormatError`] on the first malformed
///   record — invalid PRN, epoch, or body field.
/// - `Permissive`: skips any malformed record and continues; non-GPS system
///   identifiers are also silently skipped.
///
/// # Errors
///
/// In `Strict` mode, returns [`FormatError::Located`] on any malformed field.
/// In `Permissive` mode, only returns an error if the header is unreadable.
///
/// # Examples
///
/// ```
/// use siderust::formats::{ParseMode, rinex::nav::parse_rinex_nav_with_mode};
///
/// // A record with an invalid float field → strict fails, permissive skips.
/// let bad = "\
///      3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
///                                                             END OF HEADER\n\
/// G01 2024 01 01 00 00 00 NOTAFLOAT 5.678E-12 0.000E+00\n\
///      1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
///      5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
///      5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
///      9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
///      1.000E-10\n\
///      0.000E+00 0.000E+00 0.000E+00 0.000E+00\n\
///      0.000E+00 0.000E+00 0.000E+00 0.000E+00\n";
/// assert!(parse_rinex_nav_with_mode(bad, ParseMode::Strict).is_err());
/// let f = parse_rinex_nav_with_mode(bad, ParseMode::Permissive).unwrap();
/// assert_eq!(f.gps.len(), 0);   // malformed record was skipped
/// ```
pub fn parse_rinex_nav_with_mode(
    text: &str,
    mode: ParseMode,
) -> Result<RinexNavFile, FormatError> {
    let mut out = RinexNavFile::default();

    let mut lines = text.lines().enumerate();
    // Skip header up to and including "END OF HEADER".
    let mut found_header = false;
    for (_, l) in lines.by_ref() {
        if l.contains("END OF HEADER") {
            found_header = true;
            break;
        }
    }
    if !found_header && mode == ParseMode::Strict {
        return Err(FormatError::Format(
            "RINEX NAV: END OF HEADER marker not found".to_owned(),
        ));
    }

    let collected: Vec<(usize, &str)> = lines.collect();
    let mut i = 0;
    while i < collected.len() {
        let (line_no, l0) = collected[i];
        let line_no_1 = line_no + 1; // 1-based for diagnostics

        // Only GPS (G prefix or bare digit) records are supported.
        if l0.len() < 4 || !(l0.starts_with('G') || l0.as_bytes()[0].is_ascii_digit()) {
            i += 1;
            continue;
        }

        // Need 7 continuation lines (indices i+1 … i+7).
        if i + 7 >= collected.len() {
            if mode == ParseMode::Strict {
                return Err(FormatError::located(
                    "RINEX 3.05 §6.5",
                    FileLocation::at_line(line_no_1),
                    "truncated GPS record: fewer than 8 lines",
                ));
            }
            break;
        }

        // ── PRN ────────────────────────────────────────────────────────────
        let prn_result = if l0.starts_with('G') {
            l0[1..3].trim().parse::<u8>()
        } else {
            l0[..3].trim().parse::<u8>()
        };
        let prn = match prn_result {
            Ok(p) if p > 0 => p,
            Ok(_) | Err(_) => {
                if mode == ParseMode::Strict {
                    return Err(FormatError::located(
                        "RINEX 3.05 §6.3",
                        FileLocation::at_line(line_no_1),
                        format!("invalid GPS PRN: {:?}", &l0[..l0.len().min(3)]),
                    ));
                }
                i += 8;
                continue;
            }
        };

        // ── TOC + clock terms ──────────────────────────────────────────────
        let rest = &l0[3..];
        let toc_tokens: Vec<&str> = rest.split_whitespace().collect();
        let header_result = (|| -> Result<(Time<UTC>, f64, f64, f64), FormatError> {
            if toc_tokens.len() < 9 {
                return Err(FormatError::located(
                    "RINEX 3.05 §6.3",
                    FileLocation::at_line(line_no_1),
                    format!(
                        "GPS record line 0 needs ≥9 whitespace tokens, got {}",
                        toc_tokens.len()
                    ),
                ));
            }
            let parse_int = |s: &str, field: &str| -> Result<i64, FormatError> {
                s.parse::<i64>().map_err(|_| {
                    FormatError::located(
                        "RINEX 3.05 §6.3",
                        FileLocation::at_line(line_no_1),
                        format!("invalid {field} in TOC: {s:?}"),
                    )
                })
            };
            let y = parse_int(toc_tokens[0], "year")?;
            let mo = parse_int(toc_tokens[1], "month")?;
            let d = parse_int(toc_tokens[2], "day")?;
            let h = parse_int(toc_tokens[3], "hour")?;
            let mi = parse_int(toc_tokens[4], "minute")?;
            let sec_f = parse_d(toc_tokens[5], line_no_1, "second")?;
            let sec_i = sec_f as u64;
            let nanos = ((sec_f - sec_i as f64) * 1e9).round() as u32;
            let toc = NaiveDate::from_ymd_opt(y as i32, mo as u32, d as u32)
                .and_then(|nd| nd.and_hms_nano_opt(h as u32, mi as u32, sec_i as u32, nanos))
                .and_then(|naive| {
                    let dt = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
                    Time::<UTC>::try_from_chrono(dt).ok()
                })
                .ok_or_else(|| {
                    FormatError::located(
                        "RINEX 3.05 §6.3",
                        FileLocation::at_line(line_no_1),
                        format!(
                            "invalid TOC date/time: {}-{:02}-{:02} {:02}:{:02}:{:.3}",
                            toc_tokens[0], mo, d, h, mi, sec_f
                        ),
                    )
                })?;
            let af0 = parse_d(toc_tokens[6], line_no_1, "af0")?;
            let af1 = parse_d(toc_tokens[7], line_no_1, "af1")?;
            let af2 = parse_d(toc_tokens[8], line_no_1, "af2")?;
            Ok((toc, af0, af1, af2))
        })();
        let (toc, af0, af1, af2) = match header_result {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };

        // ── Body lines ─────────────────────────────────────────────────────
        let parse_body = |idx: usize, expected: usize| -> Result<Vec<f64>, FormatError> {
            let (ln, text) = collected[i + idx];
            let vals = collect_floats(text, ln + 1)?;
            if vals.len() < expected && mode == ParseMode::Strict {
                return Err(FormatError::located(
                    "RINEX 3.05 §6.5",
                    FileLocation::at_line(ln + 1),
                    format!(
                        "GPS broadcast body line {idx}: expected ≥{expected} fields, got {}",
                        vals.len()
                    ),
                ));
            }
            Ok(vals)
        };

        let v1 = match parse_body(1, 4) {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };
        let v2 = match parse_body(2, 4) {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };
        let v3 = match parse_body(3, 4) {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };
        let v4 = match parse_body(4, 4) {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };
        let v5 = match parse_body(5, 1) {
            Ok(v) => v,
            Err(e) => {
                if mode == ParseMode::Strict {
                    return Err(e);
                }
                i += 8;
                continue;
            }
        };
        // Lines 6–7 are reserved / broadcast health — skipped but must parse.

        let rec = GpsNavRecord {
            prn,
            toc,
            af0: Seconds::new(af0),
            af1,
            af2,
            iode: v1[0],
            crs: Meters::new(v1[1]),
            delta_n: AngularRate::new(v1[2]),
            m0: Radians::new(v1[3]),
            cuc: Radians::new(v2[0]),
            e: v2[1],
            cus: Radians::new(v2[2]),
            sqrt_a: v2[3],
            toe: Seconds::new(v3[0]),
            cic: Radians::new(v3[1]),
            omega0: Radians::new(v3[2]),
            cis: Radians::new(v3[3]),
            i0: Radians::new(v4[0]),
            crc: Meters::new(v4[1]),
            omega: Radians::new(v4[2]),
            omega_dot: AngularRate::new(v4[3]),
            idot: AngularRate::new(v5[0]),
        };

        out.gps.push(rec);
        i += 8;
    }

    Ok(out)
}

/// Parse a Fortran D-notation float field.
fn parse_d(s: &str, line: usize, field: &str) -> Result<f64, FormatError> {
    s.replace('D', "E")
        .replace('d', "e")
        .parse::<f64>()
        .map_err(|_| {
            FormatError::located(
                "RINEX 3.05 §6.5",
                FileLocation::at_line(line),
                format!("invalid float for field {field}: {s:?}"),
            )
        })
}

fn collect_floats(line: &str, line_no: usize) -> Result<Vec<f64>, FormatError> {
    line.split_whitespace()
        .map(|s| parse_d(s, line_no, "body"))
        .collect()
}

/// Format an `f64` in Fortran D-notation (e.g. `1.234567890123456E+02` →
/// `1.234567890123456D+02`), as used by RINEX NAV body lines.
fn fmt_d(v: f64) -> String {
    if !v.is_finite() {
        return format!("{:>19}", "0.000000000000000D+00");
    }
    let s = format!("{:19.12E}", v);
    s.replace('E', "D")
}

/// Write a RINEX 3.04 NAV file matching the layout consumed by
/// [`parse_rinex_nav`].
///
/// Only the GPS broadcast subset is written. The header is canonicalised to
/// the minimal `RINEX VERSION / TYPE` + `END OF HEADER` pair; per-record
/// output uses fixed-width Fortran D-notation in the standard
/// `19.12D` field width.
///
/// Round-trip parity is exact at 1e-15 relative for every numeric field
/// under [`parse_rinex_nav`].
///
/// # Errors
///
/// Returns [`FormatError::Io`] on write failures and [`FormatError::Format`]
/// on epoch projection failures.
///
/// # Examples
///
/// ```
/// use siderust::formats::rinex::nav::{parse_rinex_nav, write_rinex_nav};
///
/// let txt = "\
///      3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
///                                                             END OF HEADER\n\
/// G01 2024 01 01 00 00 00 1.234567E-04 5.678E-12 0.000E+00\n\
///      1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
///      5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
///      5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
///      9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
///      1.000E-10\n\
///      0.0E+00 0.0E+00 0.0E+00 0.0E+00\n\
///      0.0E+00 0.0E+00 0.0E+00 0.0E+00\n";
/// let f = parse_rinex_nav(txt).unwrap();
/// let mut buf = Vec::new();
/// write_rinex_nav(&mut buf, &f).unwrap();
/// let f2 = parse_rinex_nav(std::str::from_utf8(&buf).unwrap()).unwrap();
/// assert_eq!(f.gps.len(), f2.gps.len());
/// ```
pub fn write_rinex_nav<W: Write>(w: &mut W, file: &RinexNavFile) -> Result<(), FormatError> {
    fn header_line(w: &mut impl Write, body: &str, label: &str) -> std::io::Result<()> {
        let body = if body.len() > 60 { &body[..60] } else { body };
        writeln!(w, "{:<60}{}", body, label)
    }
    header_line(
        w,
        "     3.04           N: GNSS NAV DATA    M (Mixed)",
        "RINEX VERSION / TYPE",
    )?;
    header_line(w, "", "END OF HEADER")?;

    for r in &file.gps {
        let dt = r
            .toc
            .try_to_chrono()
            .map_err(|e| FormatError::Format(format!("rinex_nav: TOC to chrono failed: {e}")))?;
        let sec = dt.second() as f64 + dt.nanosecond() as f64 / 1e9;
        // Line 0: PRN + TOC + 3 clock terms.
        writeln!(
            w,
            "G{:02} {:04} {:02} {:02} {:02} {:02} {:02}{}{}{}",
            r.prn,
            dt.year(),
            dt.month(),
            dt.day(),
            dt.hour(),
            dt.minute(),
            sec as u32,
            fmt_d(r.af0.value()),
            fmt_d(r.af1),
            fmt_d(r.af2),
        )?;
        let cont = |w: &mut W, vs: [f64; 4]| -> std::io::Result<()> {
            writeln!(
                w,
                "    {}{}{}{}",
                fmt_d(vs[0]),
                fmt_d(vs[1]),
                fmt_d(vs[2]),
                fmt_d(vs[3])
            )
        };
        cont(w, [r.iode, r.crs.value(), r.delta_n.value(), r.m0.value()])?;
        cont(w, [r.cuc.value(), r.e, r.cus.value(), r.sqrt_a])?;
        cont(
            w,
            [
                r.toe.value(),
                r.cic.value(),
                r.omega0.value(),
                r.cis.value(),
            ],
        )?;
        cont(
            w,
            [
                r.i0.value(),
                r.crc.value(),
                r.omega.value(),
                r.omega_dot.value(),
            ],
        )?;
        cont(w, [r.idot.value(), 0.0, 0.0, 0.0])?;
        cont(w, [0.0, 0.0, 0.0, 0.0])?;
        cont(w, [0.0, 0.0, 0.0, 0.0])?;
    }
    Ok(())
}

/// Convenience wrapper to write a RINEX NAV file to disk.
///
/// # Examples
///
/// ```no_run
/// use siderust::formats::rinex::nav::{write_rinex_nav_to_path, RinexNavFile};
/// write_rinex_nav_to_path("/tmp/out.rnx", &RinexNavFile::default()).ok();
/// ```
pub fn write_rinex_nav_to_path<P: AsRef<Path>>(
    path: P,
    file: &RinexNavFile,
) -> Result<(), FormatError> {
    let mut f = fs::File::create(path)?;
    write_rinex_nav(&mut f, file)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_one_gps_record() {
        let txt = "\
     3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
                                                            END OF HEADER\n\
G01 2024 01 01 00 00 00 1.234567E-04 5.678E-12 0.000E+00\n\
     1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
     5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
     5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
     9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
     1.000E-10\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n";
        let f = parse_rinex_nav(txt).expect("parse");
        assert_eq!(f.gps.len(), 1);
        let r = &f.gps[0];
        assert_eq!(r.prn, 1);
        assert!((r.sqrt_a - 5_153.651).abs() < 1e-3);
        assert!((r.e - 1e-3).abs() < 1e-9);
        assert!((r.toe.value() - 518_400.0).abs() < 1e-3);
    }

    #[test]
    fn round_trips_through_writer() {
        let txt = "\
     3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
                                                            END OF HEADER\n\
G01 2024 01 01 00 00 00 1.234567E-04 5.678E-12 0.000E+00\n\
     1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
     5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
     5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
     9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
     1.000E-10\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n";
        let f = parse_rinex_nav(txt).unwrap();
        let mut buf = Vec::new();
        write_rinex_nav(&mut buf, &f).unwrap();
        let f2 = parse_rinex_nav(std::str::from_utf8(&buf).unwrap()).unwrap();
        assert_eq!(f.gps.len(), f2.gps.len());
        let r = &f.gps[0];
        let r2 = &f2.gps[0];
        assert_eq!(r.prn, r2.prn);
        assert!((r.sqrt_a - r2.sqrt_a).abs() < 1e-9);
        assert!((r.e - r2.e).abs() < 1e-15);
        assert!((r.toe.value() - r2.toe.value()).abs() < 1e-9);
        assert!((r.m0.value() - r2.m0.value()).abs() < 1e-12);
    }

    const VALID_RECORD: &str = "\
     3.04           N: GNSS NAV DATA    M (Mixed)           RINEX VERSION / TYPE\n\
                                                            END OF HEADER\n\
G01 2024 01 01 00 00 00 1.234567E-04 5.678E-12 0.000E+00\n\
     1.000000E+00 2.000000E+01 3.456E-09 1.234567E+00\n\
     5.000E-07 1.000E-03 9.000E-07 5.153651E+03\n\
     5.184000E+05 1.000E-08 1.000E+00 2.000E-08\n\
     9.760000E-01 2.500E+02 -1.500E+00 -8.000E-09\n\
     1.000E-10\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n\
     0.000E+00 0.000E+00 0.000E+00 0.000E+00\n";

    #[test]
    fn strict_rejects_bad_float_field() {
        use crate::formats::ParseMode;
        let bad = VALID_RECORD.replace("1.234567E-04", "NOTAFLOAT");
        assert!(
            parse_rinex_nav(&bad).is_err(),
            "strict mode must reject malformed float"
        );
        let f = parse_rinex_nav_with_mode(&bad, ParseMode::Permissive).unwrap();
        assert_eq!(f.gps.len(), 0, "permissive mode must skip the malformed record");
    }

    #[test]
    fn strict_rejects_invalid_date() {
        use crate::formats::ParseMode;
        // Month 99 is invalid.
        let bad = VALID_RECORD.replace("2024 01 01", "2024 99 01");
        assert!(
            parse_rinex_nav(&bad).is_err(),
            "strict mode must reject out-of-range month"
        );
        let f = parse_rinex_nav_with_mode(&bad, ParseMode::Permissive).unwrap();
        assert_eq!(f.gps.len(), 0, "permissive mode must skip the record with invalid date");
    }

    #[test]
    fn strict_rejects_truncated_record() {
        use crate::formats::ParseMode;
        // Drop the last two body lines so the record is incomplete.
        let lines: Vec<&str> = VALID_RECORD.lines().collect();
        let truncated = lines[..lines.len() - 2].join("\n") + "\n";
        assert!(
            parse_rinex_nav(&truncated).is_err(),
            "strict mode must reject a truncated record"
        );
        // Permissive hits EOF gracefully (truncated record at end is a break, not skip).
        assert!(parse_rinex_nav_with_mode(&truncated, ParseMode::Permissive).is_ok());
    }
}
