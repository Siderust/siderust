//! # SP3 precise orbit reader and writer
//!
//! ## Scientific scope
//!
//! SP3 files distribute precise GNSS orbit and clock products as regularly
//! sampled Cartesian state tables. This module supports the subset required
//! by the POD workspace to ingest reference trajectories and emit
//! comparable precise-orbit products.
//!
//! The implementation is intentionally format-focused. It does not
//! interpolate trajectories, reconcile multi-constellation metadata
//! differences, or estimate clocks beyond the values explicitly present in
//! the file.
//!
//! ## Technical scope
//!
//! The public APIs are `read_sp3`, `write_sp3`, and the data containers
//! `Sp3Record`, `Sp3Epoch`, and `Sp3Position`. They preserve the row-major,
//! epoch-sampled structure of the SP3 standard so service and QC code can
//! decide how to consume the product.
//!
//! Any orbit fitting, frame conversion, or residual analysis is delegated
//! to other crates.
//!
//! ## References
//!
//! - International GNSS Service. (2020). SP3-c / SP3-d Orbit Format
//!   Specification.
//! - Montenbruck, O., Steigenberger, P., & Khachikyan, R. (2017). GNSS
//!   satellite geometry and ephemeris products. GPS Solutions, 21, 101-111.
use affn::cartesian;
use affn::centers::{AffineCenter, ReferenceCenter};
use affn::frames::GCRS;
use chrono::{DateTime, NaiveDate, Utc as ChronoUtc};
use qtty::time::Microseconds;
use qtty::unit::Kilometer;

/// Geocentric center marker for SP3 positions.
#[derive(Debug, Copy, Clone)]
pub struct EarthCenter;

impl ReferenceCenter for EarthCenter {
    type Params = ();
    fn center_name() -> &'static str {
        "Geocentric"
    }
}

impl AffineCenter for EarthCenter {}

/// Local type alias for SP3 positions: geocentric, GCRS-framed, km.
type Position<F = GCRS, U = Kilometer> = cartesian::Position<EarthCenter, F, U>;
use std::io::{BufRead, BufReader, Read, Write};
use tempoch::{Time, UTC};
use thiserror::Error;

/// SP3 parse / write errors.
#[derive(Debug, Error)]
pub enum Sp3Error {
    /// I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// Header malformed.
    #[error("malformed SP3 header at line {line}: {message}")]
    Header {
        /// 1-based line number.
        line: usize,
        /// Diagnostic.
        message: String,
    },
    /// Record line malformed.
    #[error("malformed SP3 record at line {line}: {message}")]
    Record {
        /// 1-based line number.
        line: usize,
        /// Diagnostic.
        message: String,
    },
}

/// One satellite position record at one epoch (km, microsecond clock).
#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Position {
    /// Satellite identifier (3 chars, e.g. `G01`).
    pub sat_id: String,
    /// Position in GCRS frame (km).
    pub position: Position<GCRS, Kilometer>,
    /// Clock bias. 999999.999999 µs means unavailable.
    pub clock: Microseconds,
}

/// One epoch with all satellite records.
#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Epoch {
    /// UTC epoch.
    pub time: Time<UTC>,
    /// Position records.
    pub positions: Vec<Sp3Position>,
}

/// Full SP3 record.
#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Record {
    /// Verbatim header lines (preserved for round-trip).
    pub header: Vec<String>,
    /// Epochs.
    pub epochs: Vec<Sp3Epoch>,
}

fn civil_to_time_utc(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
    line: usize,
) -> Result<Time<UTC>, Sp3Error> {
    let second_int = second as u32;
    let nanos = ((second - second_int as f64) * 1e9).round() as u32;
    let naive = NaiveDate::from_ymd_opt(year, month, day)
        .and_then(|d| d.and_hms_nano_opt(hour, minute, second_int, nanos))
        .ok_or_else(|| Sp3Error::Record {
            line,
            message: "invalid epoch date/time".into(),
        })?;
    let dt: DateTime<ChronoUtc> = DateTime::from_naive_utc_and_offset(naive, ChronoUtc);
    Time::<UTC>::try_from_chrono(dt).map_err(|e| Sp3Error::Record {
        line,
        message: format!("epoch UTC conversion: {e}"),
    })
}

/// Parse an SP3 file from a reader.
pub fn read_sp3<R: Read>(r: R) -> Result<Sp3Record, Sp3Error> {
    let mut buf = BufReader::new(r);
    let mut header = Vec::new();
    let mut epochs = Vec::new();
    let mut line_no = 0usize;
    let mut in_header = true;
    let mut current: Option<Sp3Epoch> = None;

    let mut s = String::new();
    loop {
        s.clear();
        let n = buf.read_line(&mut s)?;
        if n == 0 {
            break;
        }
        line_no += 1;
        let trimmed_eol = s.trim_end_matches(['\r', '\n']).to_string();

        if trimmed_eol == "EOF" {
            break;
        }

        if in_header {
            // Header runs until the first '*' epoch line.
            if trimmed_eol.starts_with('*') {
                in_header = false;
            } else {
                header.push(trimmed_eol);
                continue;
            }
        }

        if trimmed_eol.starts_with('*') {
            // Epoch line: "* YYYY MM DD HH MM SS.SSSSSSSS"
            if let Some(e) = current.take() {
                epochs.push(e);
            }
            let fields: Vec<&str> = trimmed_eol
                .strip_prefix('*')
                .unwrap_or(trimmed_eol.as_str())
                .split_whitespace()
                .collect();
            if fields.len() < 6 {
                return Err(Sp3Error::Record {
                    line: line_no,
                    message: format!("epoch needs 6 fields, got {}", fields.len()),
                });
            }
            current = Some(Sp3Epoch {
                time: civil_to_time_utc(
                    parse_field(fields[0], line_no, "year")?,
                    parse_field(fields[1], line_no, "month")?,
                    parse_field(fields[2], line_no, "day")?,
                    parse_field(fields[3], line_no, "hour")?,
                    parse_field(fields[4], line_no, "minute")?,
                    parse_field(fields[5], line_no, "second")?,
                    line_no,
                )?,
                positions: Vec::new(),
            });
        } else if trimmed_eol.starts_with('P') {
            let epoch = current.as_mut().ok_or_else(|| Sp3Error::Record {
                line: line_no,
                message: "P record before any epoch".into(),
            })?;
            // Columns are fixed-width in SP3, but split_whitespace is
            // robust enough for the typical fixtures we ingest.
            // Format: "PXXX  X.X         Y.Y         Z.Z         CLOCK"
            let id = trimmed_eol
                .get(1..4)
                .ok_or_else(|| Sp3Error::Record {
                    line: line_no,
                    message: "P record too short for sat id".into(),
                })?
                .trim()
                .to_string();
            let rest: Vec<&str> = trimmed_eol[4..].split_whitespace().collect();
            if rest.len() < 4 {
                return Err(Sp3Error::Record {
                    line: line_no,
                    message: format!("P record needs 4 numeric fields, got {}", rest.len()),
                });
            }
            epoch.positions.push(Sp3Position {
                sat_id: id,
                position: Position::<GCRS, Kilometer>::new(
                    parse_field::<f64>(rest[0], line_no, "x")?,
                    parse_field::<f64>(rest[1], line_no, "y")?,
                    parse_field::<f64>(rest[2], line_no, "z")?,
                ),
                clock: Microseconds::new(parse_field(rest[3], line_no, "clock")?),
            });
        } else if trimmed_eol.starts_with('V')
            || trimmed_eol.starts_with("EP")
            || trimmed_eol.starts_with("EV")
        {
            // Velocity (V) and exponent (EP/EV) records are not yet
            // surfaced; ignore them silently to keep the parser
            // permissive on real-world SP3 files.
            continue;
        }
    }
    if let Some(e) = current {
        epochs.push(e);
    }
    Ok(Sp3Record { header, epochs })
}

fn parse_field<T>(raw: &str, line: usize, what: &str) -> Result<T, Sp3Error>
where
    T: std::str::FromStr,
    <T as std::str::FromStr>::Err: std::fmt::Display,
{
    raw.parse::<T>().map_err(|e| Sp3Error::Record {
        line,
        message: format!("bad {what}: {e}"),
    })
}

/// Write an SP3 record. The header lines are emitted verbatim, then each
/// epoch as `* YYYY MM DD HH MM SS.SS...` followed by `PXXX X Y Z CLK` lines.
pub fn write_sp3<W: Write>(w: &mut W, rec: &Sp3Record) -> Result<(), Sp3Error> {
    use chrono::{Datelike, Timelike};
    for h in &rec.header {
        writeln!(w, "{h}")?;
    }
    for e in &rec.epochs {
        let dt = e.time.try_to_chrono().map_err(|err| Sp3Error::Record {
            line: 0,
            message: format!("epoch to chrono: {err}"),
        })?;
        let second = dt.second() as f64 + dt.nanosecond() as f64 / 1e9;
        writeln!(
            w,
            "*  {:4} {:>2} {:>2} {:>2} {:>2} {:11.8}",
            dt.year(),
            dt.month(),
            dt.day(),
            dt.hour(),
            dt.minute(),
            second
        )?;
        for p in &e.positions {
            writeln!(
                w,
                "P{:<3} {:14.6} {:14.6} {:14.6} {:14.6}",
                p.sat_id,
                p.position.x().value(),
                p.position.y().value(),
                p.position.z().value(),
                p.clock.value()
            )?;
        }
    }
    writeln!(w, "EOF")?;
    Ok(())
}

/// Iterator-based streaming SP3 reader.
///
/// Yields one [`Sp3Epoch`] at a time from a buffered reader without
/// materialising the entire file in memory. Use this when ingesting
/// hour-of-day SP3 files that can exceed 100 MB.
///
/// The reader-borrowing iterator is a thin layer over [`read_sp3`] and
/// makes the same parsing decisions; in particular, the SP3 header is
/// pre-consumed and exposed via [`Sp3Stream::header`].
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sp3::Sp3Stream;
///
/// const SAMPLE: &str = "\
/// #dP2024  1  1  0  0  0.00000000       1 ORBIT IGS20 HLM  IGS\n\
/// ##  2295 518400.00000000   900.00000000 60310 0.0000000000000\n\
/// +    1   G01                                                      \n\
/// ++         5                                                       \n\
/// %c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n\
/// %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n\
/// %f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n\
/// %f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n\
/// %i    0    0    0    0      0      0      0      0         0\n\
/// %i    0    0    0    0      0      0      0      0         0\n\
/// /* COMMENT\n\
/// *  2024  1  1  0  0  0.00000000\n\
/// PG01  15000.000000  20000.000000  -5000.000000      0.000123\n\
/// EOF\n";
/// let mut s = Sp3Stream::new(SAMPLE.as_bytes()).unwrap();
/// assert!(!s.header().is_empty());
/// let mut epochs = 0;
/// while let Some(ep) = s.next_epoch().unwrap() {
///     assert!(!ep.positions.is_empty());
///     epochs += 1;
/// }
/// assert_eq!(epochs, 1);
/// ```
pub struct Sp3Stream<R: Read> {
    inner: BufReader<R>,
    header: Vec<String>,
    pending_epoch_line: Option<(usize, String)>,
    line_no: usize,
    finished: bool,
}

impl<R: Read> Sp3Stream<R> {
    /// Create a new streaming reader by consuming the SP3 header up to the
    /// first epoch line.
    pub fn new(r: R) -> Result<Self, Sp3Error> {
        let mut inner = BufReader::new(r);
        let mut header = Vec::new();
        let mut line_no = 0usize;
        let mut s = String::new();
        let pending = loop {
            s.clear();
            let n = inner.read_line(&mut s)?;
            if n == 0 {
                break None;
            }
            line_no += 1;
            let trimmed = s.trim_end_matches(['\r', '\n']).to_string();
            if trimmed == "EOF" {
                break None;
            }
            if trimmed.starts_with('*') {
                break Some((line_no, trimmed));
            }
            header.push(trimmed);
        };
        Ok(Self {
            inner,
            header,
            pending_epoch_line: pending,
            line_no,
            finished: false,
        })
    }

    /// Header lines collected before the first epoch.
    pub fn header(&self) -> &[String] {
        &self.header
    }

    /// Yield the next [`Sp3Epoch`], or `None` when the stream is exhausted.
    pub fn next_epoch(&mut self) -> Result<Option<Sp3Epoch>, Sp3Error> {
        if self.finished {
            return Ok(None);
        }
        let (epoch_line_no, epoch_line) = match self.pending_epoch_line.take() {
            Some(v) => v,
            None => {
                self.finished = true;
                return Ok(None);
            }
        };
        let mut epoch = parse_epoch_line(&epoch_line, epoch_line_no)?;
        let mut s = String::new();
        loop {
            s.clear();
            let n = self.inner.read_line(&mut s)?;
            if n == 0 {
                self.finished = true;
                break;
            }
            self.line_no += 1;
            let trimmed = s.trim_end_matches(['\r', '\n']).to_string();
            if trimmed == "EOF" {
                self.finished = true;
                break;
            }
            if trimmed.starts_with('*') {
                self.pending_epoch_line = Some((self.line_no, trimmed));
                break;
            }
            if trimmed.starts_with('P') {
                epoch.positions.push(parse_p_line(&trimmed, self.line_no)?);
            }
            // V/EP/EV records ignored as in batch reader.
        }
        Ok(Some(epoch))
    }
}

fn parse_epoch_line(line: &str, line_no: usize) -> Result<Sp3Epoch, Sp3Error> {
    let fields: Vec<&str> = line
        .strip_prefix('*')
        .unwrap_or(line)
        .split_whitespace()
        .collect();
    if fields.len() < 6 {
        return Err(Sp3Error::Record {
            line: line_no,
            message: format!("epoch needs 6 fields, got {}", fields.len()),
        });
    }
    Ok(Sp3Epoch {
        time: civil_to_time_utc(
            parse_field(fields[0], line_no, "year")?,
            parse_field(fields[1], line_no, "month")?,
            parse_field(fields[2], line_no, "day")?,
            parse_field(fields[3], line_no, "hour")?,
            parse_field(fields[4], line_no, "minute")?,
            parse_field(fields[5], line_no, "second")?,
            line_no,
        )?,
        positions: Vec::new(),
    })
}

fn parse_p_line(line: &str, line_no: usize) -> Result<Sp3Position, Sp3Error> {
    let id = line
        .get(1..4)
        .ok_or_else(|| Sp3Error::Record {
            line: line_no,
            message: "P record too short for sat id".into(),
        })?
        .trim()
        .to_string();
    let rest: Vec<&str> = line[4..].split_whitespace().collect();
    if rest.len() < 4 {
        return Err(Sp3Error::Record {
            line: line_no,
            message: format!("P record needs 4 numeric fields, got {}", rest.len()),
        });
    }
    Ok(Sp3Position {
        sat_id: id,
        position: Position::<GCRS, Kilometer>::new(
            parse_field::<f64>(rest[0], line_no, "x")?,
            parse_field::<f64>(rest[1], line_no, "y")?,
            parse_field::<f64>(rest[2], line_no, "z")?,
        ),
        clock: Microseconds::new(parse_field(rest[3], line_no, "clock")?),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
#dP2024  1  1  0  0  0.00000000      96 ORBIT IGS20 HLM  IGS\n\
##  2295 518400.00000000   900.00000000 60310 0.0000000000000\n\
+    2   G01G02                                                     \n\
++         5   5                                                    \n\
%c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n\
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n\
%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n\
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n\
%i    0    0    0    0      0      0      0      0         0\n\
%i    0    0    0    0      0      0      0      0         0\n\
/* COMMENT LINE\n\
*  2024  1  1  0  0  0.00000000\n\
PG01   1000.000000   2000.000000   3000.000000      0.000123\n\
PG02  -1500.000000   1700.000000  -2500.000000      0.000456\n\
*  2024  1  1  0 15  0.00000000\n\
PG01   1100.000000   2100.000000   3100.000000      0.000124\n\
PG02  -1400.000000   1800.000000  -2400.000000      0.000457\n\
EOF\n";

    #[test]
    fn parse_sample() {
        let rec = read_sp3(SAMPLE.as_bytes()).expect("parse sample");
        assert_eq!(rec.epochs.len(), 2);
        assert_eq!(rec.epochs[0].positions.len(), 2);
        assert_eq!(rec.epochs[0].positions[0].sat_id, "G01");
        assert!((rec.epochs[0].positions[0].position.x().value() - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn round_trip() {
        let rec = read_sp3(SAMPLE.as_bytes()).expect("parse");
        let mut buf = Vec::new();
        write_sp3(&mut buf, &rec).expect("write");
        let rec2 = read_sp3(buf.as_slice()).expect("re-parse");
        assert_eq!(rec.epochs.len(), rec2.epochs.len());
        assert_eq!(rec.header.len(), rec2.header.len());
        for (e1, e2) in rec.epochs.iter().zip(rec2.epochs.iter()) {
            // Epoch times may differ by ~14 µs due to TT_MINUS_TAI imprecision in
            // the try_to_chrono backward path; tolerate up to 1 ms.
            let t1 = e1.time.raw_seconds_pair().0.value();
            let t2 = e2.time.raw_seconds_pair().0.value();
            assert!((t1 - t2).abs() < 1e-3, "epoch time mismatch: {t1} vs {t2}");
            assert_eq!(e1.positions, e2.positions);
        }
    }
}
