//! # IGS Long File Name (LFN) helpers
//!
//! ## Scientific scope
//!
//! The International GNSS Service defines a Long File Name convention that
//! encodes analysis centre, solution type, time span, and sample rate into
//! a single deterministic filename. Correct naming is essential for IGS
//! product archives, downstream automation, and multi-agency product
//! comparisons.
//!
//! This module provides helpers for SP3, RINEX observation, and clock file
//! naming. It does **not** validate the scientific content of the files.
//!
//! ## Technical scope
//!
//! The IGS LFN format (v2) is:
//!
//! ```text
//! {AAA}{V}{NNNNNN}_{YYYY}{DOY}{HH}{MM}_{LEN}_{SMP}_{FTY}.{EXT}
//! ```
//!
//! - `AAA`    – 3-char analysis centre (e.g. `IGS`, `GFZ`, `EUR`).
//! - `V`      – 1-char version (`0`–`9` or `a`–`z`).
//! - `NNNNNN` – 6-char campaign/solution token (e.g. `OPSFIN`, `OPSRAP`).
//! - `YYYY`   – 4-digit year.
//! - `DOY`    – 3-digit day-of-year (001–366).
//! - `HH`     – 2-digit UTC hour (00–23).
//! - `MM`     – 2-digit UTC minute (00–59).
//! - `LEN`    – 3-char time span, e.g. `01D`, `24H`, `15M`.
//! - `SMP`    – 3-char sample rate, e.g. `15M`, `30S`, `05S`.
//! - `FTY`    – 3-char file type, e.g. `ORB`, `CLK`, `ATT`.
//! - `EXT`    – file extension without leading dot, e.g. `SP3`, `CLK`, `rnx`.
//!
//! ## References
//!
//! - International GNSS Service. (2020). IGS File and Product Naming
//!   Conventions. Version 2.
//! - International GNSS Service. (2020). SP3-c / SP3-d Orbit Format
//!   Specification.

use chrono::{Datelike, NaiveDate};
use std::path::PathBuf;

// ─── SP3 type/version markers ───────────────────────────────────────────────

/// SP3 file content type: position-only or position-and-velocity.
///
/// Determines the IGS LFN file-type tag and the SP3 header type character.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::Sp3FileType;
///
/// assert_eq!(Sp3FileType::PositionOnly.lfn_file_type_tag(), "ORB");
/// assert_eq!(Sp3FileType::PositionVelocity.lfn_file_type_tag(), "OBV");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sp3FileType {
    /// File contains position records only (`P` records; SP3 header type `P`).
    PositionOnly,
    /// File contains position and velocity records (`P`+`V`; SP3 header type `V`).
    PositionVelocity,
}

impl Sp3FileType {
    /// Returns the 3-character IGS LFN file-type tag for this content type.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::naming::Sp3FileType;
    ///
    /// assert_eq!(Sp3FileType::PositionOnly.lfn_file_type_tag(), "ORB");
    /// ```
    pub fn lfn_file_type_tag(self) -> &'static str {
        match self {
            Self::PositionOnly => "ORB",
            Self::PositionVelocity => "OBV",
        }
    }
}

/// SP3 format version marker.
///
/// Identifies which revision of the SP3 standard a file conforms to.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::Sp3Version;
///
/// assert_eq!(Sp3Version::D.lfn_version_char(), 'd');
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sp3Version {
    /// SP3-a — legacy format; header `#a`.
    A,
    /// SP3-c — current standard; header `#c`.
    C,
    /// SP3-d — latest standard; header `#d`.
    D,
}

impl Sp3Version {
    /// Returns the version character used in the IGS LFN `V` field.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::naming::Sp3Version;
    ///
    /// assert_eq!(Sp3Version::C.lfn_version_char(), 'c');
    /// ```
    pub fn lfn_version_char(self) -> char {
        match self {
            Self::A => 'a',
            Self::C => 'c',
            Self::D => 'd',
        }
    }
}

// ─── Sample-rate encoding ────────────────────────────────────────────────────

/// Encode a sample interval (seconds) as a 3-character IGS LFN sample-rate token.
///
/// | Range | Format | Example |
/// |---|---|---|
/// | `< 60 s` | `NNS` | `30` → `"30S"` |
/// | `60 – 3 599 s` | `NNM` | `900` → `"15M"` |
/// | `≥ 3 600 s` | `NNH` | `3600` → `"01H"` |
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::sample_rate_to_lfn;
///
/// assert_eq!(sample_rate_to_lfn(30),   "30S");
/// assert_eq!(sample_rate_to_lfn(900),  "15M");
/// assert_eq!(sample_rate_to_lfn(3600), "01H");
/// ```
pub fn sample_rate_to_lfn(sample_rate_sec: u32) -> String {
    if sample_rate_sec < 60 {
        format!("{:02}S", sample_rate_sec)
    } else if sample_rate_sec < 3600 {
        format!("{:02}M", sample_rate_sec / 60)
    } else {
        format!("{:02}H", sample_rate_sec / 3600)
    }
}

// ─── High-level SP3 naming ───────────────────────────────────────────────────

/// Build an IGS LFN SP3 orbit product filename from high-level parameters.
///
/// Encodes the analysis-centre code (`satellite`), calendar `date`,
/// `file_type`, `sample_rate_sec`, and `version` into the IGS Long File
/// Name convention. The campaign token defaults to `"OPSFIN"` (operational
/// final), span to `"01D"` (one day), and start time to `00:00` UTC.
///
/// Use [`sp3_lfn`] for full control over every LFN field.
///
/// # Examples
///
/// ```
/// use chrono::NaiveDate;
/// use siderust::pod::product::naming::{Sp3FileType, Sp3Version, sp3_filename};
///
/// let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
/// let p = sp3_filename("IGS", date, Sp3FileType::PositionOnly, 900, Sp3Version::D);
/// assert_eq!(
///     p.file_name().unwrap().to_str().unwrap(),
///     "IGSdOPSFIN_20240010000_01D_15M_ORB.SP3",
/// );
/// ```
pub fn sp3_filename(
    satellite: &str,
    date: NaiveDate,
    file_type: Sp3FileType,
    sample_rate_sec: u32,
    version: Sp3Version,
) -> PathBuf {
    let year = date.year();
    let doy = date.ordinal();
    let v = version.lfn_version_char();
    let smp = sample_rate_to_lfn(sample_rate_sec);
    let fty = file_type.lfn_file_type_tag();
    PathBuf::from(format!(
        "{satellite}{v}OPSFIN_{year:04}{doy:03}0000_01D_{smp}_{fty}.SP3"
    ))
}

/// Build an IGS short-name SP3 filename (pre-LFN convention).
///
/// The traditional IGS short name encodes the GPS week and day-of-week:
/// `{ac}{gpsweek:04}{dow}.sp3`. The GPS epoch is 1980-01-06 (Sunday);
/// `dow` counts from Sunday (0) to Saturday (6).
///
/// Use [`sp3_filename`] for the current IGS Long File Name convention.
///
/// # Examples
///
/// ```
/// use chrono::NaiveDate;
/// use siderust::pod::product::naming::{Sp3Version, sp3_short_filename};
///
/// // 2024-01-01 falls on GPS week 2295, day-of-week 1 (Monday).
/// let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
/// let p = sp3_short_filename("igs", date, Sp3Version::D);
/// assert_eq!(p.to_str().unwrap(), "igs22951.sp3");
/// ```
pub fn sp3_short_filename(ac: &str, date: NaiveDate, _version: Sp3Version) -> PathBuf {
    let gps_epoch = NaiveDate::from_ymd_opt(1980, 1, 6).expect("GPS epoch is valid");
    let days = (date - gps_epoch).num_days();
    let gps_week = days / 7;
    let dow = days % 7;
    PathBuf::from(format!("{ac}{gps_week:04}{dow}.sp3"))
}

// ─── Full-control LFN builder ────────────────────────────────────────────────

/// Build a generic IGS Long File Name with full control over every field.
///
/// Prefer the typed wrappers ([`sp3_filename`], [`sp3_lfn`], [`clk_filename`],
/// [`rnx_obs_filename`]) for routine use.
///
/// # Arguments
///
/// * `ac`        – Analysis centre, exactly 3 ASCII chars (e.g. `"IGS"`).
/// * `version`   – Single version digit or letter (e.g. `'0'`).
/// * `campaign`  – Campaign / solution token, exactly 6 ASCII chars (e.g. `"OPSFIN"`).
/// * `year`      – 4-digit year.
/// * `doy`       – 1-based day-of-year (1–366).
/// * `hour`      – UTC hour (0–23).
/// * `minute`    – UTC minute (0–59).
/// * `span`      – Time span string, exactly 3 chars (e.g. `"01D"`).
/// * `sample`    – Sample rate string, exactly 3 chars (e.g. `"15M"`).
/// * `file_type` – File-type tag, exactly 3 chars (e.g. `"ORB"`).
/// * `ext`       – File extension without leading dot (e.g. `"SP3"`).
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::igs_lfn;
///
/// let p = igs_lfn("IGS", '0', "OPSFIN", 2024, 1, 0, 0, "01D", "15M", "ORB", "SP3");
/// assert_eq!(p.file_name().unwrap().to_str().unwrap(),
///            "IGS0OPSFIN_20240010000_01D_15M_ORB.SP3");
/// ```
#[allow(clippy::too_many_arguments)]
pub fn igs_lfn(
    ac: &str,
    version: char,
    campaign: &str,
    year: i32,
    doy: u32,
    hour: u32,
    minute: u32,
    span: &str,
    sample: &str,
    file_type: &str,
    ext: &str,
) -> PathBuf {
    PathBuf::from(format!(
        "{ac}{version}{campaign}_{year:04}{doy:03}{hour:02}{minute:02}_{span}_{sample}_{file_type}.{ext}",
    ))
}

/// Build an IGS SP3 orbit product filename with full LFN field control.
///
/// Convenience wrapper around [`igs_lfn`] that fixes `file_type` to `"ORB"`
/// and `ext` to `"SP3"`. Use [`sp3_filename`] for the simpler high-level API.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::sp3_lfn;
///
/// let p = sp3_lfn("IGS", '0', "OPSFIN", 2024, 1, 0, 0, "01D", "15M");
/// assert_eq!(
///     p.file_name().unwrap().to_str().unwrap(),
///     "IGS0OPSFIN_20240010000_01D_15M_ORB.SP3",
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn sp3_lfn(
    ac: &str,
    version: char,
    campaign: &str,
    year: i32,
    doy: u32,
    hour: u32,
    minute: u32,
    span: &str,
    sample: &str,
) -> PathBuf {
    igs_lfn(
        ac, version, campaign, year, doy, hour, minute, span, sample, "ORB", "SP3",
    )
}

/// Build an IGS clock product filename.
///
/// Convenience wrapper that fixes `file_type` to `"CLK"` and `ext` to `"CLK"`.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::clk_filename;
///
/// let p = clk_filename("GFZ", '0', "OPSRAP", 2024, 1, 0, 0, "01D", "30S");
/// assert_eq!(
///     p.file_name().unwrap().to_str().unwrap(),
///     "GFZ0OPSRAP_20240010000_01D_30S_CLK.CLK",
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn clk_filename(
    ac: &str,
    version: char,
    campaign: &str,
    year: i32,
    doy: u32,
    hour: u32,
    minute: u32,
    span: &str,
    sample: &str,
) -> PathBuf {
    igs_lfn(
        ac, version, campaign, year, doy, hour, minute, span, sample, "CLK", "CLK",
    )
}

/// Build an IGS RINEX observation filename.
///
/// Convenience wrapper that fixes `file_type` to `"MO"` (mixed observation)
/// and `ext` to `"rnx"` (RINEX 3+ lowercase convention).
///
/// # Examples
///
/// ```
/// use siderust::pod::product::naming::rnx_obs_filename;
///
/// let p = rnx_obs_filename("EUR", '0', "MGXFIN", 2024, 1, 0, 0, "01D", "30S");
/// assert_eq!(
///     p.file_name().unwrap().to_str().unwrap(),
///     "EUR0MGXFIN_20240010000_01D_30S_MO.rnx",
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn rnx_obs_filename(
    ac: &str,
    version: char,
    campaign: &str,
    year: i32,
    doy: u32,
    hour: u32,
    minute: u32,
    span: &str,
    sample: &str,
) -> PathBuf {
    igs_lfn(
        ac, version, campaign, year, doy, hour, minute, span, sample, "MO", "rnx",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn igs_lfn_format_matches_spec() {
        let p = igs_lfn(
            "IGS", '0', "OPSFIN", 2024, 1, 0, 0, "01D", "15M", "ORB", "SP3",
        );
        assert_eq!(
            p.to_str().unwrap(),
            "IGS0OPSFIN_20240010000_01D_15M_ORB.SP3"
        );
    }

    #[test]
    fn sp3_lfn_matches_igs_lfn() {
        let a = sp3_lfn("GFZ", '0', "OPSRAP", 2024, 32, 12, 0, "01D", "05M");
        let b = igs_lfn(
            "GFZ", '0', "OPSRAP", 2024, 32, 12, 0, "01D", "05M", "ORB", "SP3",
        );
        assert_eq!(a, b);
    }

    #[test]
    fn sp3_filename_high_level() {
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let p = sp3_filename("IGS", date, Sp3FileType::PositionOnly, 900, Sp3Version::D);
        assert_eq!(
            p.to_str().unwrap(),
            "IGSdOPSFIN_20240010000_01D_15M_ORB.SP3"
        );
    }

    #[test]
    fn sp3_filename_velocity_type() {
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let p = sp3_filename(
            "GFZ",
            date,
            Sp3FileType::PositionVelocity,
            30,
            Sp3Version::C,
        );
        let s = p.to_str().unwrap();
        assert!(s.contains("OBV"), "expected OBV file-type tag, got {s}");
        assert!(s.contains("30S"), "expected 30S sample rate, got {s}");
    }

    #[test]
    fn sp3_short_filename_gps_week() {
        // 2024-01-01 = GPS week 2295, DOW 1 (Monday).
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let p = sp3_short_filename("igs", date, Sp3Version::D);
        assert_eq!(p.to_str().unwrap(), "igs22951.sp3");
    }

    #[test]
    fn sp3_short_filename_gps_epoch() {
        // GPS epoch itself = GPS week 0, DOW 0 (Sunday).
        let date = NaiveDate::from_ymd_opt(1980, 1, 6).unwrap();
        let p = sp3_short_filename("igs", date, Sp3Version::A);
        assert_eq!(p.to_str().unwrap(), "igs00000.sp3");
    }

    #[test]
    fn clk_filename_correct_extension() {
        let p = clk_filename("GFZ", '0', "OPSRAP", 2024, 1, 0, 0, "01D", "30S");
        let s = p.to_str().unwrap();
        assert!(s.ends_with(".CLK"), "expected .CLK, got {s}");
    }

    #[test]
    fn rnx_obs_filename_correct_extension() {
        let p = rnx_obs_filename("EUR", '0', "MGXFIN", 2024, 1, 0, 0, "01D", "30S");
        let s = p.to_str().unwrap();
        assert!(s.ends_with(".rnx"), "expected .rnx, got {s}");
    }

    #[test]
    fn doy_zero_padded() {
        let p = sp3_lfn("IGS", '0', "OPSFIN", 2024, 7, 0, 0, "01D", "15M");
        assert!(
            p.to_str().unwrap().contains("2024007"),
            "DOY must be zero-padded to 3 digits"
        );
    }

    #[test]
    fn sample_rate_encoding() {
        assert_eq!(sample_rate_to_lfn(5), "05S");
        assert_eq!(sample_rate_to_lfn(30), "30S");
        assert_eq!(sample_rate_to_lfn(60), "01M");
        assert_eq!(sample_rate_to_lfn(900), "15M");
        assert_eq!(sample_rate_to_lfn(3600), "01H");
    }
}
