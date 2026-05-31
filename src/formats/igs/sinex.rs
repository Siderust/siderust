//! # SINEX station-coordinate reader
//!
//! Parses the Solution INdependent EXchange (SINEX) format to extract
//! station coordinate and velocity estimates from geodetic solutions.
//!
//! Supported blocks: `+SITE/ID`, `+SOLUTION/ESTIMATE`.
//! Estimate types: `STAX`, `STAY`, `STAZ` (position), `VELX`, `VELY`, `VELZ`
//! (velocity).
//!
//! ## References
//!
//! - IERS/IGS SINEX Format Description Version 2.10.

use super::{FileLocation, FormatError, ParseMode};
use affn::cartesian;
use affn::centers::{AffineCenter, ReferenceCenter};
use affn::frames::ITRF;
use qtty::unit::Meter;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

/// Earth geocenter marker for SINEX ITRF coordinates.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sinex::GeoCenterItrf;
/// use affn::centers::ReferenceCenter;
/// assert_eq!(GeoCenterItrf::center_name(), "Geocentric ITRF");
/// ```
#[derive(Debug, Copy, Clone)]
pub struct GeoCenterItrf;

impl ReferenceCenter for GeoCenterItrf {
    type Params = ();
    fn center_name() -> &'static str {
        "Geocentric ITRF"
    }
}

impl AffineCenter for GeoCenterItrf {}

/// SINEX ITRF Cartesian position (meters).
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sinex::SinexPosition;
/// let _pos = SinexPosition::new(4_641_949.0_f64, 1_393_045.0_f64, 4_133_287.0_f64);
/// ```
pub type SinexPosition = cartesian::Position<GeoCenterItrf, ITRF, Meter>;

/// Velocity unit: meters per year.
///
/// Used for VELX/VELY/VELZ estimates in SINEX files.
pub type MeterPerYear = qtty::Per<qtty::unit::Meter, qtty::unit::Year>;

/// A station solution extracted from a SINEX file.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sinex::{StationCoordinate, SinexPosition};
/// let pos = SinexPosition::new(1.0_f64, 2.0_f64, 3.0_f64);
/// let _ = StationCoordinate {
///     code: "MATE".to_string(),
///     point_code: "A".to_string(),
///     solution_id: "1".to_string(),
///     ref_epoch_mjd: 51544.0,
///     position: pos,
///     velocity_m_yr: None,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct StationCoordinate {
    /// 4-character station code (DOMES notation).
    pub code: String,
    /// Point code (e.g. `"A"`).
    pub point_code: String,
    /// Solution ID (e.g. `"1"`).
    pub solution_id: String,
    /// Reference epoch (MJD, f64).
    pub ref_epoch_mjd: f64,
    /// ITRF Cartesian position (meters).
    pub position: SinexPosition,
    /// Velocity vector components in meters/year, if present.
    pub velocity_m_yr: Option<[f64; 3]>,
}

/// Parsed SINEX solution.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sinex::SinexSolution;
/// let sol = SinexSolution::default();
/// assert!(sol.stations.is_empty());
/// ```
#[derive(Debug, Default)]
pub struct SinexSolution {
    /// File creation agency.
    pub agency: String,
    /// Solution description.
    pub description: String,
    /// Station coordinate/velocity records.
    pub stations: Vec<StationCoordinate>,
}

/// Read a SINEX file from a byte source.
///
/// Parses `+SITE/ID`, `+SOLUTION/EPOCHS`, and `+SOLUTION/ESTIMATE` blocks.
/// Estimates for `STAX`/`STAY`/`STAZ` become position components;
/// `VELX`/`VELY`/`VELZ` become velocity components.
///
/// # Errors
///
/// Returns [`FormatError::Format`] for fatal parse errors in `Strict` mode,
/// or [`FormatError::Io`] for underlying I/O failures.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::sinex::read_sinex;
/// use siderust::formats::ParseMode;
///
/// let data = b"%=SNX 2.10 IGS 24:001:00000 IGS 00:000:00000 23:365:86370 P 00000 0 S\n\
///              +SITE/ID\n\
///              *Code Pt Domes____ T _Station Description__ _Longitude_ _Latitude__ _Height_\n\
///              ABCD  A 10101M001 P Some Station            010 00 00.0 050 00 00.0  100.000\n\
///              -SITE/ID\n\
///              +SOLUTION/ESTIMATE\n\
///              *Index Type__ Code Pt Soln _Ref_Epoch__ Unit S ___Estimated_Value___ __Std_Dev__\n\
///              1 STAX   ABCD  A    1 00:001:00000 m    2  1.000000000000000E+06  1.000E-04\n\
///              2 STAY   ABCD  A    1 00:001:00000 m    2 -2.000000000000000E+06  1.000E-04\n\
///              3 STAZ   ABCD  A    1 00:001:00000 m    2  3.000000000000000E+06  1.000E-04\n\
///              -SOLUTION/ESTIMATE\n\
///              %ENDSNX\n";
///
/// let sol = read_sinex(&data[..], ParseMode::Permissive).unwrap();
/// assert_eq!(sol.stations.len(), 1);
/// assert_eq!(sol.stations[0].code, "ABCD");
/// ```
pub fn read_sinex<R: Read>(reader: R, mode: ParseMode) -> Result<SinexSolution, FormatError> {
    let mut solution = SinexSolution::default();
    let buf = BufReader::new(reader);

    // Accumulate SOLUTION/ESTIMATE rows keyed by station code.
    #[derive(Default)]
    struct EstRow {
        xyz: [Option<f64>; 3],
        vel: [Option<f64>; 3],
        ref_epoch_mjd: f64,
        pt: String,
        soln: String,
    }
    let mut estimates: HashMap<String, EstRow> = HashMap::new();

    #[derive(PartialEq)]
    enum Block {
        None,
        SiteId,
        Estimate,
    }
    let mut block = Block::None;

    for (lineno, result) in buf.lines().enumerate() {
        let line = result.map_err(FormatError::Io)?;
        let lineno = lineno + 1;

        // SINEX header
        if line.starts_with("%=SNX") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                solution.agency = parts[2].to_string();
            }
            continue;
        }
        if line.starts_with("%ENDSNX") {
            break;
        }
        // Comment lines inside blocks
        if line.starts_with('*') || line.starts_with('%') {
            continue;
        }
        // Block start/end markers
        if let Some(bname) = line.strip_prefix('+') {
            match bname.trim() {
                "SITE/ID" => block = Block::SiteId,
                "SOLUTION/ESTIMATE" => block = Block::Estimate,
                _ => block = Block::None,
            }
            continue;
        }
        if line.starts_with('-') {
            block = Block::None;
            continue;
        }

        match block {
            Block::SiteId => {
                // We only record that the site exists; ignore detailed site info.
                // The station code is in the first token.
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.is_empty() {
                    continue;
                }
                // Just note the code exists; coordinate data comes from ESTIMATE block.
                let _ = parts[0];
            }
            Block::Estimate => {
                // Format: index type code pt soln epoch unit constraint value stddev
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() < 9 {
                    if mode == ParseMode::Strict {
                        return Err(FormatError::located(
                            "SINEX 2.10 §4.2",
                            FileLocation::at_line(lineno),
                            format!("ESTIMATE line has {} fields, need ≥ 9", parts.len()),
                        ));
                    }
                    continue;
                }
                let est_type = parts[1];
                let code = parts[2].to_string();
                let pt = parts[3].to_string();
                let soln = parts[4].to_string();
                let epoch_raw = parts[5];
                let value: f64 = match parts[8].parse() {
                    Ok(v) => v,
                    Err(_) => {
                        if mode == ParseMode::Strict {
                            return Err(FormatError::located(
                                "SINEX 2.10 §4.2",
                                FileLocation::at_line(lineno),
                                format!("cannot parse value: {}", parts[8]),
                            ));
                        }
                        continue;
                    }
                };

                let ref_epoch_mjd = sinex_epoch_to_mjd(epoch_raw);
                let row = estimates.entry(code.clone()).or_default();
                row.pt = pt;
                row.soln = soln;
                row.ref_epoch_mjd = ref_epoch_mjd;
                match est_type {
                    "STAX" => row.xyz[0] = Some(value),
                    "STAY" => row.xyz[1] = Some(value),
                    "STAZ" => row.xyz[2] = Some(value),
                    "VELX" => row.vel[0] = Some(value),
                    "VELY" => row.vel[1] = Some(value),
                    "VELZ" => row.vel[2] = Some(value),
                    _ => {}
                }
            }
            Block::None => {}
        }
    }

    // Build station records from accumulated estimates.
    for (code, row) in estimates {
        let (x, y, z) = match (row.xyz[0], row.xyz[1], row.xyz[2]) {
            (Some(x), Some(y), Some(z)) => (x, y, z),
            _ => continue,
        };
        let velocity_m_yr = match (row.vel[0], row.vel[1], row.vel[2]) {
            (Some(vx), Some(vy), Some(vz)) => Some([vx, vy, vz]),
            _ => None,
        };
        solution.stations.push(StationCoordinate {
            code,
            point_code: row.pt,
            solution_id: row.soln,
            ref_epoch_mjd: row.ref_epoch_mjd,
            position: SinexPosition::new(x, y, z),
            velocity_m_yr,
        });
    }

    Ok(solution)
}

/// Convert a SINEX epoch string `"YY:DDD:SSSSS"` to MJD.
fn sinex_epoch_to_mjd(s: &str) -> f64 {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return 0.0;
    }
    let yy: i32 = parts[0].parse().unwrap_or(0);
    let doy: u32 = parts[1].parse().unwrap_or(1);
    let sod: f64 = parts[2].parse().unwrap_or(0.0);
    let yyyy = if yy == 0 {
        2000
    } else if yy >= 80 {
        1900 + yy
    } else {
        2000 + yy
    };
    let epoch = match chrono::NaiveDate::from_yo_opt(yyyy, doy.max(1)) {
        Some(d) => d,
        None => return 0.0,
    };
    let mjd_ref = chrono::NaiveDate::from_ymd_opt(1858, 11, 17).expect("valid MJD epoch");
    let days = (epoch - mjd_ref).num_days() as f64;
    days + sod / 86400.0
}
