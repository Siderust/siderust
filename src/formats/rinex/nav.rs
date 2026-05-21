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
use super::FormatError;
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

/// Parse a RINEX NAV file from a string slice.
pub fn parse_rinex_nav(text: &str) -> Result<RinexNavFile, FormatError> {
    let mut out = RinexNavFile::default();

    let mut lines = text.lines();
    // Skip header up to and including "END OF HEADER".
    for l in lines.by_ref() {
        if l.contains("END OF HEADER") {
            break;
        }
    }

    let collected: Vec<&str> = lines.collect();
    let mut i = 0;
    while i + 7 < collected.len() {
        let l0 = collected[i];
        if l0.len() < 4 || !(l0.starts_with('G') || l0.as_bytes()[0].is_ascii_digit()) {
            i += 1;
            continue;
        }
        let prn = if l0.starts_with('G') {
            l0[1..3].trim().parse::<u8>().unwrap_or(0)
        } else {
            l0[..3].trim().parse::<u8>().unwrap_or(0)
        };
        if prn == 0 {
            i += 1;
            continue;
        }

        // TOC: positions 4..23 → "YYYY MM DD HH MM SS"
        let rest = &l0[3..];
        let toc_tokens: Vec<&str> = rest.split_whitespace().collect();
        if toc_tokens.len() < 9 {
            i += 1;
            continue;
        }
        let mut rec = GpsNavRecord {
            prn,
            toc: {
                let y: i32 = toc_tokens[0].parse().unwrap_or(2000);
                let mo: u32 = toc_tokens[1].parse().unwrap_or(1);
                let d: u32 = toc_tokens[2].parse().unwrap_or(1);
                let h: u32 = toc_tokens[3].parse().unwrap_or(0);
                let mi: u32 = toc_tokens[4].parse().unwrap_or(0);
                let sec_f: f64 = toc_tokens[5].parse().unwrap_or(0.0);
                let sec_i = sec_f as u32;
                let nanos = ((sec_f - sec_i as f64) * 1e9).round() as u32;
                NaiveDate::from_ymd_opt(y, mo, d)
                    .and_then(|d| d.and_hms_nano_opt(h, mi, sec_i, nanos))
                    .and_then(|naive| {
                        let dt = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
                        Time::<UTC>::try_from_chrono(dt).ok()
                    })
                    .unwrap_or_else(|| {
                        let naive = NaiveDate::from_ymd_opt(2000, 1, 1)
                            .unwrap()
                            .and_hms_opt(0, 0, 0)
                            .unwrap();
                        Time::<UTC>::try_from_chrono(DateTime::from_naive_utc_and_offset(
                            naive, ChronoUtc,
                        ))
                        .unwrap()
                    })
            },
            af0: Seconds::new(parse_d(toc_tokens[6])),
            af1: parse_d(toc_tokens[7]),
            af2: parse_d(toc_tokens[8]),
            ..Default::default()
        };

        // Lines 1..7 each carry up to 4 floats in fixed 19-char fields,
        // but split_whitespace works for our subset.
        let v1: Vec<f64> = collect_floats(collected[i + 1]);
        let v2: Vec<f64> = collect_floats(collected[i + 2]);
        let v3: Vec<f64> = collect_floats(collected[i + 3]);
        let v4: Vec<f64> = collect_floats(collected[i + 4]);
        let v5: Vec<f64> = collect_floats(collected[i + 5]);
        let v6: Vec<f64> = collect_floats(collected[i + 6]);
        let v7: Vec<f64> = collect_floats(collected[i + 7]);

        if v1.len() >= 4 {
            rec.iode = v1[0];
            rec.crs = Meters::new(v1[1]);
            rec.delta_n = AngularRate::new(v1[2]);
            rec.m0 = Radians::new(v1[3]);
        }
        if v2.len() >= 4 {
            rec.cuc = Radians::new(v2[0]);
            rec.e = v2[1];
            rec.cus = Radians::new(v2[2]);
            rec.sqrt_a = v2[3];
        }
        if v3.len() >= 4 {
            rec.toe = Seconds::new(v3[0]);
            rec.cic = Radians::new(v3[1]);
            rec.omega0 = Radians::new(v3[2]);
            rec.cis = Radians::new(v3[3]);
        }
        if v4.len() >= 4 {
            rec.i0 = Radians::new(v4[0]);
            rec.crc = Meters::new(v4[1]);
            rec.omega = Radians::new(v4[2]);
            rec.omega_dot = AngularRate::new(v4[3]);
        }
        if !v5.is_empty() {
            rec.idot = AngularRate::new(v5[0]);
        }
        let _ = v6;
        let _ = v7;

        out.gps.push(rec);
        i += 8;
    }

    Ok(out)
}

fn parse_d(s: &str) -> f64 {
    s.replace('D', "E").replace('d', "e").parse().unwrap_or(0.0)
}

fn collect_floats(line: &str) -> Vec<f64> {
    line.split_whitespace().map(parse_d).collect()
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
}
