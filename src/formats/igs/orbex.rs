//! # ORBEX orbit/clock/attitude product reader
//!
//! ORBEX is a flexible, text-based format for distributing precise satellite
//! orbit, clock, and attitude data. This module parses the `#ORB`, `#CLK`,
//! and `#ATT` record types from the `+EPHEMERIS/DATA` block.
//!
//! ## References
//!
//! - IGS ORBEX Format Description (draft, 2020).

use super::{FileLocation, FormatError, ParseMode};
use std::io::{BufRead, BufReader, Read};

/// Single orbit record from an ORBEX `#ORB` block.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::orbex::OrbexOrbitEntry;
/// let e = OrbexOrbitEntry {
///     sat: "G01".to_string(),
///     epoch_mjd: 60000.25,
///     pos_m: [1e7, 2e6, 3e6],
///     vel_m_s: None,
/// };
/// assert_eq!(e.sat, "G01");
/// ```
#[derive(Debug, Clone)]
pub struct OrbexOrbitEntry {
    /// Satellite identifier (e.g. `"G01"`).
    pub sat: String,
    /// Epoch as modified Julian date.
    pub epoch_mjd: f64,
    /// Position components in meters.
    pub pos_m: [f64; 3],
    /// Velocity components in m/s, if present.
    pub vel_m_s: Option<[f64; 3]>,
}

/// Single clock record from an ORBEX `#CLK` block.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::orbex::OrbexClockEntry;
/// let e = OrbexClockEntry { sat: "G01".to_string(), epoch_mjd: 60000.25, bias_s: 1e-7 };
/// assert!((e.bias_s - 1e-7).abs() < 1e-15);
/// ```
#[derive(Debug, Clone)]
pub struct OrbexClockEntry {
    /// Satellite identifier.
    pub sat: String,
    /// Epoch as modified Julian date.
    pub epoch_mjd: f64,
    /// Clock bias in seconds.
    pub bias_s: f64,
}

/// Single attitude record from an ORBEX `#ATT` block.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::orbex::OrbexAttitudeEntry;
/// let e = OrbexAttitudeEntry {
///     sat: "G01".to_string(),
///     epoch_mjd: 60000.25,
///     quaternion: [1.0, 0.0, 0.0, 0.0],
/// };
/// assert_eq!(e.quaternion[0], 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct OrbexAttitudeEntry {
    /// Satellite identifier.
    pub sat: String,
    /// Epoch as modified Julian date.
    pub epoch_mjd: f64,
    /// Quaternion \[w, x, y, z\].
    pub quaternion: [f64; 4],
}

/// Parsed ORBEX product.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::orbex::OrbexProduct;
/// let prod = OrbexProduct::default();
/// assert!(prod.orbits.is_empty());
/// ```
#[derive(Debug, Default)]
pub struct OrbexProduct {
    /// File version string.
    pub version: String,
    /// List of orbit entries.
    pub orbits: Vec<OrbexOrbitEntry>,
    /// List of clock entries.
    pub clocks: Vec<OrbexClockEntry>,
    /// List of attitude entries.
    pub attitudes: Vec<OrbexAttitudeEntry>,
}

/// Read an ORBEX file from a byte source.
///
/// Parses `#ORB`, `#CLK`, and `#ATT` record types. Unknown record types are
/// skipped in both `Strict` and `Permissive` modes (they are a defined
/// extension mechanism in the ORBEX spec).
///
/// # Errors
///
/// Returns [`FormatError::Format`] on parse failure in `Strict` mode,
/// or [`FormatError::Io`] for I/O failures.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::orbex::read_orbex;
/// use siderust::formats::ParseMode;
///
/// let lines = [
///     "%=ORBEX  0.09", "%%", "+EPHEMERIS/DATA", "#ORB",
///     "## 2024  1 21600.000000000",
///     " G01  1.0E+07  2.0E+06  3.0E+06",
///     "-EPHEMERIS/DATA", "%ENDORBEX",
/// ];
/// let data = lines.join("\n");
/// let prod = read_orbex(data.as_bytes(), ParseMode::Permissive).unwrap();
/// assert_eq!(prod.orbits.len(), 1);
/// ```
pub fn read_orbex<R: Read>(reader: R, mode: ParseMode) -> Result<OrbexProduct, FormatError> {
    let mut product = OrbexProduct::default();
    let buf = BufReader::new(reader);

    #[derive(PartialEq, Clone, Copy)]
    enum Section {
        None,
        Orb,
        Clk,
        Att,
    }

    let mut section = Section::None;
    let mut current_epoch_mjd: f64 = 0.0;

    for (lineno, result) in buf.lines().enumerate() {
        let line = result.map_err(FormatError::Io)?;
        let lineno = lineno + 1;

        if line.starts_with("%=ORBEX") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                product.version = parts[1].to_string();
            }
            continue;
        }
        if line.starts_with("%ENDORBEX") || line.starts_with("%EOF") {
            break;
        }
        // Skip comment/description block lines
        if line.starts_with('%') {
            continue;
        }
        if line.starts_with('+') || line.starts_with('-') {
            continue;
        }

        // Epoch header: ## YYYY DDD SSSSS.SSSSSSS (must check before generic '#')
        if line.starts_with("##") {
            let rest = line.trim_start_matches('#').trim();
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.len() >= 3 {
                let year: i32 = parts[0].parse().unwrap_or(2000);
                let doy: u32 = parts[1].parse().unwrap_or(1);
                let sod: f64 = parts[2].parse().unwrap_or(0.0);
                current_epoch_mjd = orbex_epoch_to_mjd(year, doy, sod);
            }
            continue;
        }

        // Section type markers
        if line.starts_with("#ORB") {
            section = Section::Orb;
            continue;
        }
        if line.starts_with("#CLK") {
            section = Section::Clk;
            continue;
        }
        if line.starts_with("#ATT") {
            section = Section::Att;
            continue;
        }
        // Other # markers: switch to none
        if line.starts_with('#') {
            section = Section::None;
            continue;
        }

        // Satellite data lines start with a space then sat id
        if !line.starts_with(' ') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }

        match section {
            Section::Orb => {
                if parts.len() < 4 {
                    if mode == ParseMode::Strict {
                        return Err(FormatError::located(
                            "ORBEX §3.3",
                            FileLocation::at_line(lineno),
                            format!("ORB record has {} fields, need ≥ 4", parts.len()),
                        ));
                    }
                    continue;
                }
                let sat = parts[0].to_string();
                let x = parse_f64(parts[1], lineno, "pos_x", mode)?;
                let y = parse_f64(parts[2], lineno, "pos_y", mode)?;
                let z = parse_f64(parts[3], lineno, "pos_z", mode)?;
                let vel_m_s = if parts.len() >= 7 {
                    let vx = parse_f64(parts[4], lineno, "vel_x", mode)?;
                    let vy = parse_f64(parts[5], lineno, "vel_y", mode)?;
                    let vz = parse_f64(parts[6], lineno, "vel_z", mode)?;
                    Some([vx, vy, vz])
                } else {
                    None
                };
                product.orbits.push(OrbexOrbitEntry {
                    sat,
                    epoch_mjd: current_epoch_mjd,
                    pos_m: [x, y, z],
                    vel_m_s,
                });
            }
            Section::Clk => {
                if parts.len() < 2 {
                    if mode == ParseMode::Strict {
                        return Err(FormatError::located(
                            "ORBEX §3.4",
                            FileLocation::at_line(lineno),
                            format!("CLK record has {} fields, need ≥ 2", parts.len()),
                        ));
                    }
                    continue;
                }
                let sat = parts[0].to_string();
                let bias_s = parse_f64(parts[1], lineno, "bias", mode)?;
                product.clocks.push(OrbexClockEntry {
                    sat,
                    epoch_mjd: current_epoch_mjd,
                    bias_s,
                });
            }
            Section::Att => {
                if parts.len() < 5 {
                    if mode == ParseMode::Strict {
                        return Err(FormatError::located(
                            "ORBEX §3.5",
                            FileLocation::at_line(lineno),
                            format!("ATT record has {} fields, need ≥ 5", parts.len()),
                        ));
                    }
                    continue;
                }
                let sat = parts[0].to_string();
                let w = parse_f64(parts[1], lineno, "q_w", mode)?;
                let x = parse_f64(parts[2], lineno, "q_x", mode)?;
                let y = parse_f64(parts[3], lineno, "q_y", mode)?;
                let z = parse_f64(parts[4], lineno, "q_z", mode)?;
                product.attitudes.push(OrbexAttitudeEntry {
                    sat,
                    epoch_mjd: current_epoch_mjd,
                    quaternion: [w, x, y, z],
                });
            }
            Section::None => {}
        }
    }

    Ok(product)
}

fn parse_f64(s: &str, lineno: usize, what: &str, mode: ParseMode) -> Result<f64, FormatError> {
    s.parse::<f64>().map_err(|_| {
        if mode == ParseMode::Strict {
            FormatError::located(
                "ORBEX §3",
                FileLocation::at_line(lineno),
                format!("cannot parse {what}: {s:?}"),
            )
        } else {
            FormatError::Format(format!("ORBEX line {lineno}: cannot parse {what}: {s:?}"))
        }
    })
}

fn orbex_epoch_to_mjd(year: i32, doy: u32, sod: f64) -> f64 {
    let epoch = match chrono::NaiveDate::from_yo_opt(year, doy.max(1)) {
        Some(d) => d,
        None => return 0.0,
    };
    let mjd_ref = chrono::NaiveDate::from_ymd_opt(1858, 11, 17).expect("valid MJD epoch");
    let days = (epoch - mjd_ref).num_days() as f64;
    days + sod / 86400.0
}
