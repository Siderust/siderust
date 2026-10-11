//! # CPF orbit prediction reader and writer
//!
//! ## Scientific scope
//!
//! The Consolidated Prediction Format (CPF) distributes tabulated satellite
//! orbit predictions used by SLR stations to point their telescopes and to
//! validate POD solutions. This module supports CPF versions 1 and 2 in
//! both strict and permissive parse modes.
//!
//! Position records (type `10`) are mapped to typed
//! [`CpfEphemerisEntry`] values carrying:
//!
//! - A [`tempoch::Time<UTC>`] epoch computed from the MJD + seconds-of-day.
//! - A Cartesian position in the ITRF frame expressed in kilometres as
//!   [`affn::cartesian::Position<GeocentricCenter, ITRF, Kilometer>`].
//!
//! Velocity records (type `20`), manoeuvre records, and other optional
//! sections are outside the current scope.
//!
//! The [`write_cpf`] function serialises a [`CpfFile`] back to the CPF
//! text format, enabling round-trip validation.
//!
//! ## Technical scope
//!
//! Public entry points are [`read_cpf`], [`parse_cpf`], and
//! [`parse_cpf_with_mode`]. All three return a [`CpfFile`] whose
//! tabulated ephemeris is the typed [`CpfEphemeris`] view.
//!
//! CPF v2 specification: <https://ilrs.gsfc.nasa.gov/docs/2018/cpf_2.00h-1.pdf>.
//!
//! ## References
//!
//! - International Laser Ranging Service. (2018). Consolidated Prediction
//!   Format (CPF) Specification (v2.00).
//!   <https://ilrs.gsfc.nasa.gov/docs/2018/cpf_2.0h-1.pdf>
//! - Pearlman, M. R., Degnan, J. J., & Bosworth, J. M. (2002). The
//!   International Laser Ranging Service. Advances in Space Research,
//!   30(2), 135–143.
use super::{FileLocation, FormatError, ParseMode};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar};
use affn::cartesian;
use affn::centers::{AffineCenter, ReferenceCenter};
use affn::frames::ITRF;
use alloc::format;
use alloc::string::{String, ToString};
use alloc::vec::Vec;
use qtty::unit::Kilometer;
#[cfg(feature = "std")]
use std::fs;
#[cfg(feature = "std")]
use std::io::Write;
#[cfg(feature = "std")]
use std::path::Path;
use tempoch::{Time, UTC};

// ── GeocentricCenter ──────────────────────────────────────────────────────────

/// Geocentric reference center for CPF positions (Earth's centre of mass).
///
/// CPF positions are always given relative to the geocentre in ITRF. This
/// marker type anchors the type-level center constraint for
/// [`CpfEphemerisEntry::position`].
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::GeocentricCenter;
/// use affn::centers::ReferenceCenter;
/// assert_eq!(GeocentricCenter::center_name(), "Geocentric");
/// ```
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub struct GeocentricCenter;

impl ReferenceCenter for GeocentricCenter {
    type Params = ();
    fn center_name() -> &'static str {
        "Geocentric"
    }
}

impl AffineCenter for GeocentricCenter {}

/// Typed CPF position: geocentric, ITRF frame, kilometres.
type CpfCartesian = cartesian::Position<GeocentricCenter, ITRF, Kilometer>;

// ── CpfEphemerisEntry ─────────────────────────────────────────────────────────

/// A typed CPF position record (type `10`).
///
/// Carries a UTC epoch and a geocentric ITRF position in kilometres,
/// derived from the file-native metres and MJD + SOD values.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::{parse_cpf, GeocentricCenter};
/// use affn::centers::ReferenceCenter;
///
/// let txt = "H1 CPF 2 HTS 2024 01 01 00 1\n\
///            H2 lageos1 1155 7603901 2024 1 1 0 0 0 2024 1 1 1 0 0 60 0 ITRF2014\n\
///            10 1 60310 0.0 0 6000000.0 9000000.0 7000000.0\n\
///            99\n";
/// let f = parse_cpf(txt).unwrap();
/// let e = &f.ephemeris.entries[0];
/// // Position is in km; file gives metres, so 6 000 000 m → 6 000 km
/// assert!((e.position.x().value() - 6_000.0).abs() < 1e-6);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CpfEphemerisEntry {
    /// UTC epoch (MJD midnight + seconds of day).
    pub epoch: Time<UTC>,
    /// Geocentric ITRF position in kilometres.
    pub position: CpfCartesian,
}

// ── CpfEphemeris ──────────────────────────────────────────────────────────────

/// Collection of typed CPF position records.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::{parse_cpf, CpfEphemeris};
///
/// let txt = "H1 CPF 2 CNES 2024 01 01 00 1\n\
///            H2 lageos1 1155 7603901 2024 1 1 0 0 0 2024 1 1 1 0 0 60 0 ITRF2014\n\
///            10 1 60310 0.0 0 6000000.0 9000000.0 7000000.0\n\
///            10 1 60310 60.0 0 6001200.0 8999100.0 7001500.0\n\
///            99\n";
/// let f = parse_cpf(txt).unwrap();
/// assert_eq!(f.ephemeris.entries.len(), 2);
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CpfEphemeris {
    /// Typed position entries in chronological order.
    pub entries: Vec<CpfEphemerisEntry>,
}

impl CpfEphemeris {
    /// Iterate over all ephemeris entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::ilrs::cpf::parse_cpf;
    ///
    /// let txt = "H1 CPF 2 HTS 2024 01 01 00 1\n\
    ///            H2 lageos1 1155 7603901 2024 1 1 0 0 0 2024 1 1 0 1 0 60 0 ITRF2014\n\
    ///            10 1 60310 0.0 0 7000000.0 0.0 0.0\n\
    ///            99\n";
    /// let f = parse_cpf(txt).unwrap();
    /// for e in f.ephemeris.iter() {
    ///     assert!(e.position.x().value().abs() > 1_000.0); // in km
    /// }
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = &CpfEphemerisEntry> {
        self.entries.iter()
    }
}

impl core::ops::Index<usize> for CpfEphemeris {
    type Output = CpfEphemerisEntry;

    fn index(&self, i: usize) -> &CpfEphemerisEntry {
        &self.entries[i]
    }
}

// ── CpfFile ───────────────────────────────────────────────────────────────────

/// Parsed CPF file.
///
/// The tabulated orbit lives in [`ephemeris`](CpfFile::ephemeris) as
/// typed [`CpfEphemerisEntry`] values (UTC epoch, ITRF kilometres).
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::parse_cpf;
/// let txt = "H1 CPF 2 CNES 2024 01 01 00 1\n\
///            H2 lageos1 1155 7603901 2024 1 1 0 0 0 2024 1 1 1 0 0 60 0 ITRF2014\n\
///            10 1 60310 0.0 0 6000000.0 9000000.0 7000000.0\n99\n";
/// let f = parse_cpf(txt).unwrap();
/// assert_eq!(f.version, "2");
/// assert_eq!(f.reference_frame, "ITRF2014");
/// assert_eq!(f.ephemeris.entries.len(), 1);
/// ```
#[derive(Debug, Clone, Default)]
pub struct CpfFile {
    /// Source agency (H1).
    pub source: String,
    /// CPF version (H1).
    pub version: String,
    /// Target name (H2).
    pub target_name: String,
    /// COSPAR ID (H2).
    pub cospar_id: String,
    /// SIC (H2).
    pub sic: i32,
    /// NORAD (H2).
    pub norad: String,
    /// Reference frame string from H2 (e.g. `"ITRF2014"`).
    pub reference_frame: String,
    /// Step size in seconds (H2, 0 if unspecified).
    pub time_step_s: f64,
    /// Typed ephemeris entries (km, UTC epoch).
    pub ephemeris: CpfEphemeris,
    /// Parse mode used.
    pub parse_mode: ParseMode,
}

// ── Public entry points ───────────────────────────────────────────────────────

/// Read a CPF file from disk using [`ParseMode::Strict`].
///
/// # Examples
///
/// ```no_run
/// use siderust::formats::ilrs::cpf::read_cpf;
/// let f = read_cpf("lageos1.cpf").unwrap();
/// println!("{} ephemeris entries", f.ephemeris.entries.len());
/// ```
pub fn read_cpf<P: AsRef<Path>>(path: P) -> Result<CpfFile, FormatError> {
    let text = fs::read_to_string(path.as_ref())?;
    parse_cpf_impl(&text, ParseMode::Strict, Some(path.as_ref().into()))
}

/// Parse a CPF string using [`ParseMode::Strict`].
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::parse_cpf;
/// let txt = "H1 CPF 2 HTS 2024 01 01 00 1\n\
///            H2 s 0 0 2024 1 1 0 0 0 2024 1 1 0 1 0 60 0 ITRF2014\n\
///            10 1 60310 0.0 0 7000000.0 0.0 0.0\n99\n";
/// let f = parse_cpf(txt).unwrap();
/// assert_eq!(f.version, "2");
/// ```
pub fn parse_cpf(text: &str) -> Result<CpfFile, FormatError> {
    parse_cpf_impl(text, ParseMode::Strict, None)
}

/// Parse a CPF string with an explicit [`ParseMode`].
///
/// In `Permissive` mode, malformed position records are silently skipped
/// rather than causing a hard error.
///
/// # Examples
///
/// ```
/// use siderust::formats::{ParseMode, ilrs::cpf::parse_cpf_with_mode};
/// let txt = "H1 CPF 2 HTS 2024 01 01 00 1\nH2 s 0 0\n10 1 BROKEN\n99\n";
/// let f = parse_cpf_with_mode(txt, ParseMode::Permissive).unwrap();
/// assert_eq!(f.ephemeris.entries.len(), 0); // skipped
/// ```
pub fn parse_cpf_with_mode(text: &str, mode: ParseMode) -> Result<CpfFile, FormatError> {
    parse_cpf_impl(text, mode, None)
}

/// Write a [`CpfFile`] to a writer in CPF text format.
///
/// The output is deterministic and can be re-parsed by [`parse_cpf`] to
/// recover the same records (round-trip). Positions are written as
/// metres to preserve precision.
///
/// # Examples
///
/// ```
/// use siderust::formats::ilrs::cpf::{parse_cpf, write_cpf};
///
/// let txt = "H1 CPF 2 CNES 2024 01 01 00 1\n\
///            H2 lageos1 1155 7603901 2024 1 1 0 0 0 2024 1 1 1 0 0 60 0 ITRF2014\n\
///            10 1 60310 0.000 0 6000000.000 9000000.000 7000000.000\n\
///            10 1 60310 60.000 0 6001200.000 8999100.000 7001500.000\n\
///            99\n";
/// let f = parse_cpf(txt).unwrap();
/// let mut buf = Vec::new();
/// write_cpf(&f, &mut buf).unwrap();
/// let s = String::from_utf8(buf).unwrap();
/// let f2 = parse_cpf(&s).unwrap();
/// assert_eq!(f.ephemeris.entries.len(), f2.ephemeris.entries.len());
/// assert_eq!(f.reference_frame, f2.reference_frame);
/// ```
pub fn write_cpf<W: Write>(f: &CpfFile, w: &mut W) -> Result<(), FormatError> {
    let ctx = tempoch::TimeContext::new();
    // H1: CPF <version> <source> <year> <month> <day> <hour> 1
    writeln!(w, "H1 CPF {} {} 2024 01 01 00 1", f.version, f.source).map_err(FormatError::Io)?;
    // H2: <cospar> <sic> <norad> <start> <end> <step> <leap> <frame>
    writeln!(
        w,
        "H2 {} {} {} 2024 01 01 00 00 00 2024 01 01 01 00 00 {} 0 {}",
        if f.cospar_id.is_empty() {
            &f.target_name
        } else {
            &f.cospar_id
        },
        f.sic,
        if f.norad.is_empty() { "0" } else { &f.norad },
        if f.time_step_s > 0.0 {
            format!("{}", f.time_step_s as u64)
        } else {
            "60".to_string()
        },
        f.reference_frame,
    )
    .map_err(FormatError::Io)?;
    // Position records.
    for entry in &f.ephemeris.entries {
        use chrono::Timelike;
        let civil = entry
            .epoch
            .try_to_chrono_with(&ctx)
            .map_err(|e| FormatError::Format(format!("CPF 10: UTC encoding: {e}")))?;
        let midnight = civil
            .date_naive()
            .and_time(chrono::NaiveTime::MIN)
            .and_utc();
        let mjd_int = Time::<UTC>::try_from_chrono_with(midnight, &ctx)
            .and_then(|epoch| epoch.to_with::<tempoch::MJD>(&ctx))
            .and_then(|epoch| epoch.try_raw_with(&ctx))
            .map_err(|e| FormatError::Format(format!("CPF 10: UTC date: {e}")))?
            .value() as i64;
        let sod =
            civil.num_seconds_from_midnight() as f64 + civil.nanosecond() as f64 / 1_000_000_000.0;
        // Convert km back to metres for the on-the-wire representation.
        let x_m = entry.position.x().value() * 1_000.0;
        let y_m = entry.position.y().value() * 1_000.0;
        let z_m = entry.position.z().value() * 1_000.0;
        writeln!(
            w,
            "10 1 {} {:.9} 0 {:.3} {:.3} {:.3}",
            mjd_int, sod, x_m, y_m, z_m,
        )
        .map_err(FormatError::Io)?;
    }
    writeln!(w, "99").map_err(FormatError::Io)?;
    Ok(())
}

// ── Internal parser ───────────────────────────────────────────────────────────

fn parse_cpf_impl(
    text: &str,
    mode: ParseMode,
    path: Option<std::path::PathBuf>,
) -> Result<CpfFile, FormatError> {
    let ctx = tempoch::TimeContext::new();
    let mut out = CpfFile {
        parse_mode: mode,
        ..Default::default()
    };

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
                // H1 CPF <version> <source> <year> <month> <day> <hour> <seq>
                let _fmt = tokens.next(); // "CPF"
                if let Some(v) = tokens.next() {
                    out.version = v.to_string();
                }
                if let Some(s) = tokens.next() {
                    out.source = s.to_string();
                }
            }
            "H2" => {
                // v2: H2 cospar sic norad start... end... step leap frame
                // v1: H2 cospar sic norad start_year ... frame
                // Collect all tokens; reference frame is the last one.
                let toks: Vec<&str> = tokens.collect();
                if let Some(cospar) = toks.first() {
                    out.cospar_id = cospar.to_string();
                    out.target_name = cospar.to_string();
                }
                if let Some(sic_str) = toks.get(1) {
                    out.sic = sic_str.parse().unwrap_or(0);
                }
                if let Some(norad_str) = toks.get(2) {
                    out.norad = norad_str.to_string();
                }
                // Time step: field index 15 for v2 (0-indexed), else scan.
                if let Some(step_str) = toks.get(15) {
                    out.time_step_s = step_str.parse().unwrap_or(0.0);
                }
                // Reference frame: last token.
                if let Some(frame) = toks.last() {
                    out.reference_frame = frame.to_string();
                }
            }
            "10" => {
                // 10 <dir> <mjd_int> <sod> <leap> <x_m> <y_m> <z_m>
                let rest: Vec<&str> = tokens.collect();

                let maybe_mjd = rest.get(1).and_then(|s| s.parse::<i64>().ok());
                let maybe_sod = rest.get(2).and_then(|s| s.parse::<f64>().ok());
                let maybe_x = rest.get(4).and_then(|s| s.parse::<f64>().ok());
                let maybe_y = rest.get(5).and_then(|s| s.parse::<f64>().ok());
                let maybe_z = rest.get(6).and_then(|s| s.parse::<f64>().ok());

                let (mjd_i, sod, x_m, y_m, z_m) =
                    match (maybe_mjd, maybe_sod, maybe_x, maybe_y, maybe_z) {
                        (Some(m), Some(s), Some(x), Some(y), Some(z)) => (m, s, x, y, z),
                        (None, ..) => {
                            let loc = FileLocation::new(
                                path.as_ref().map(|p| p.display().to_string()),
                                Some(line_no),
                                None,
                            );
                            let err =
                                FormatError::located("CPF v2 §4.1", loc, "record 10: missing MJD");
                            if mode == ParseMode::Strict {
                                return Err(err);
                            }
                            continue;
                        }
                        (_, None, ..) => {
                            let loc = FileLocation::new(
                                path.as_ref().map(|p| p.display().to_string()),
                                Some(line_no),
                                None,
                            );
                            let err =
                                FormatError::located("CPF v2 §4.1", loc, "record 10: missing SOD");
                            if mode == ParseMode::Strict {
                                return Err(err);
                            }
                            continue;
                        }
                        _ => {
                            let loc = FileLocation::new(
                                path.as_ref().map(|p| p.display().to_string()),
                                Some(line_no),
                                None,
                            );
                            let err = FormatError::located(
                                "CPF v2 §4.1",
                                loc,
                                "record 10: missing X/Y/Z",
                            );
                            if mode == ParseMode::Strict {
                                return Err(err);
                            }
                            continue;
                        }
                    };

                // Typed record: convert metres to km, compute epoch.
                let epoch = build_epoch_from_mjd_sod(mjd_i, sod, &ctx)?;
                out.ephemeris.entries.push(CpfEphemerisEntry {
                    epoch,
                    position: CpfCartesian::new(x_m / 1_000.0, y_m / 1_000.0, z_m / 1_000.0),
                });
            }
            "20" | "30" | "40" | "50" | "60" => { /* velocity, corrections, etc. — skipped */ }
            "99" => break, // end of file
            _ => {}        // unknown tags silently ignored
        }
    }

    Ok(out)
}

/// Build a `Time<UTC>` from an integer MJD and fractional seconds-of-day.
fn build_epoch_from_mjd_sod(
    mjd_int: i64,
    sod: f64,
    ctx: &tempoch::TimeContext,
) -> Result<Time<UTC>, FormatError> {
    // Only midnight is encoded as MJD. The time-of-day stays a civil label,
    // including chrono's extended nanosecond representation of 23:59:60.
    if !sod.is_finite() || !(0.0..86_401.0).contains(&sod) {
        return Err(FormatError::Format(format!(
            "CPF 10: invalid civil seconds-of-day {sod}"
        )));
    }
    let day_start =
        tempoch::ModifiedJulianDate::<UTC>::try_new_with(qtty::Day::new(mjd_int as f64), ctx)
            .map(|epoch| epoch.to_j2000s())
            .map_err(|e| FormatError::Format(format!("CPF 10: UTC date: {e}")))?;
    let midnight = day_start
        .try_to_chrono_with(ctx)
        .map_err(|e| FormatError::Format(format!("CPF 10: UTC date: {e}")))?;
    let whole = sod.floor() as u32;
    let second = whole.min(86_399);
    let nanos = ((sod - second as f64) * 1_000_000_000.0).round() as u32;
    let time = chrono::NaiveTime::from_num_seconds_from_midnight_opt(second, nanos)
        .ok_or_else(|| FormatError::Format("CPF 10: invalid civil time".into()))?;
    Time::<UTC>::try_from_chrono_with(midnight.date_naive().and_time(time).and_utc(), ctx)
        .map_err(|e| FormatError::Format(format!("CPF 10: UTC conversion: {e}")))
}

#[cfg(test)]
mod tests {
    use super::ParseMode;
    use super::*;
    use chrono::Timelike;

    const SAMPLE: &str = "\
H1 CPF 2 HTS 2024 01 01 00 1\n\
H2 lageos1 1155 7603901 2024 01 01 00 00 00 0 0 ITRF2014\n\
10 1 60310 0.000     0  7000000.0      0.0      0.0\n\
10 1 60310 60.000    0  7000100.0   1000.0    -50.0\n\
99\n";

    #[test]
    fn parses_minimal_cpf() {
        let f = parse_cpf(SAMPLE).expect("parse");
        assert_eq!(f.version, "2");
        assert_eq!(f.reference_frame, "ITRF2014");
        assert_eq!(f.ephemeris.entries.len(), 2);
        // 7 000 100 m → 7 000.1 km in the typed ephemeris.
        assert!(
            (f.ephemeris.entries[1].position.x().value() - 7_000.1).abs() < 1e-6,
            "X position parsed"
        );
    }

    #[test]
    fn ephemeris_entries_populated() {
        let f = parse_cpf(SAMPLE).expect("parse");
        assert_eq!(f.ephemeris.entries.len(), 2);
        let e0 = &f.ephemeris.entries[0];
        // 7 000 000 m → 7 000 km
        assert!((e0.position.x().value() - 7_000.0).abs() < 1e-6);
        assert!((e0.position.y().value()).abs() < 1e-6);
    }

    #[test]
    fn ephemeris_epoch_correct() {
        let f = parse_cpf(SAMPLE).expect("parse");
        let dt = f.ephemeris.entries[0].epoch.try_to_chrono().unwrap();
        // MJD 60310 = 2024-01-01 UTC, SOD = 0 → midnight
        assert_eq!(dt.hour(), 0);
        assert_eq!(dt.minute(), 0);
        assert_eq!(dt.second(), 0);
    }

    #[test]
    fn ephemeris_epoch_sod_offset() {
        let f = parse_cpf(SAMPLE).expect("parse");
        let dt = f.ephemeris.entries[1].epoch.try_to_chrono().unwrap();
        // SOD = 60 s → 00:01:00
        assert_eq!(dt.hour(), 0);
        assert_eq!(dt.minute(), 1);
        assert_eq!(dt.second(), 0);
    }

    #[test]
    fn geocentric_center_name() {
        use affn::centers::ReferenceCenter;
        assert_eq!(GeocentricCenter::center_name(), "Geocentric");
    }

    #[test]
    fn round_trip_write_reparse() {
        let txt = "\
H1 CPF 2 CNES 2024 01 01 00 1\n\
H2 lageos1 1155 7603901 2024 01 01 00 00 00 2024 01 01 01 00 00 60 0 ITRF2014\n\
10 1 60310 0.000 0 6000000.000 9000000.000 7000000.000\n\
10 1 60310 60.000 0 6001200.000 8999100.000 7001500.000\n\
99\n";
        let f1 = parse_cpf(txt).expect("parse original");
        let mut buf = Vec::new();
        write_cpf(&f1, &mut buf).expect("write");
        let s = String::from_utf8(buf).expect("utf8");
        let f2 = parse_cpf(&s).expect("re-parse");
        assert_eq!(
            f1.ephemeris.entries.len(),
            f2.ephemeris.entries.len(),
            "ephemeris entry count preserved"
        );
        assert_eq!(
            f1.reference_frame, f2.reference_frame,
            "reference frame preserved"
        );
        for (p1, p2) in f1.ephemeris.entries.iter().zip(f2.ephemeris.entries.iter()) {
            assert!(
                (p1.position.x().value() - p2.position.x().value()).abs() < 1e-6,
                "X position round-trip"
            );
            assert!(
                (p1.position.y().value() - p2.position.y().value()).abs() < 1e-6,
                "Y position round-trip"
            );
            assert!(
                (p1.position.z().value() - p2.position.z().value()).abs() < 1e-6,
                "Z position round-trip"
            );
        }
    }

    #[test]
    fn strict_missing_mjd_fails() {
        let bad = "H1 CPF 2 TST 2024 01 01 00 1\nH2 tgt 0 0\n10 1\n";
        let err = parse_cpf(bad).expect_err("missing MJD should fail");
        assert!(matches!(err, FormatError::Located { .. }));
    }

    #[test]
    fn permissive_missing_mjd_skips() {
        let bad = "H1 CPF 2 TST 2024 01 01 00 1\nH2 tgt 0 0\n10 1\n";
        let f = parse_cpf_with_mode(bad, ParseMode::Permissive).expect("permissive ok");
        assert_eq!(f.ephemeris.entries.len(), 0);
    }

    #[test]
    fn strict_missing_xyz_fails() {
        let bad = "H1 CPF 2 TST 2024 01 01 00 1\nH2 tgt 0 0\n10 1 60310 0.0 0\n";
        let err = parse_cpf(bad).expect_err("missing XYZ should fail");
        assert!(matches!(err, FormatError::Located { .. }));
    }

    #[test]
    fn permissive_missing_xyz_skips() {
        let bad = "H1 CPF 2 TST 2024 01 01 00 1\nH2 tgt 0 0\n10 1 60310 0.0 0\n";
        let f = parse_cpf_with_mode(bad, ParseMode::Permissive).expect("permissive ok");
        assert_eq!(f.ephemeris.entries.len(), 0);
    }

    #[test]
    fn unknown_tags_ignored() {
        let txt = "H1 CPF 2 CNES 2024 01 01 00 1\nXX garbage\nH2 sat 0 0\n99\n";
        let f = parse_cpf(txt).expect("no error");
        assert_eq!(f.source, "CNES");
    }
}
