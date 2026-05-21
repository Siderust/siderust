// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Column-positional parser for the canonical NORAD 2LE / 3LE format.

use chrono::{Duration, NaiveDate, TimeZone, Utc};
use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;
use tempoch::{Time, UTC};

use super::TleError;
use crate::formats::tle::{Classification, InternationalDesignator, SatelliteNumber, TLE};

/// Validate the canonical TLE checksum (digits + minus signs, mod 10).
///
/// The checksum character is the last column (column 69) of the line.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::validate_tle_checksum;
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// validate_tle_checksum(l1).unwrap();
/// ```
pub fn validate_tle_checksum(line: &str) -> Result<(), TleError> {
    let line_no = if line.starts_with('1') {
        1u8
    } else if line.starts_with('2') {
        2u8
    } else {
        0u8
    };
    if line.len() != 69 {
        return Err(TleError::BadLength {
            line: line_no,
            got: line.len(),
        });
    }
    let bytes = line.as_bytes();
    let stated_byte = bytes[68];
    if !stated_byte.is_ascii_digit() {
        return Err(TleError::BadChecksum {
            line: line_no,
            stated: 0,
            computed: compute_checksum(&line[..68]),
        });
    }
    let stated = stated_byte - b'0';
    let computed = compute_checksum(&line[..68]);
    if stated != computed {
        return Err(TleError::BadChecksum {
            line: line_no,
            stated,
            computed,
        });
    }
    Ok(())
}

/// Compute the canonical TLE checksum digit for the first 68 columns.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::compute_tle_checksum;
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  292";
/// assert_eq!(compute_tle_checksum(l1), 7);
/// ```
pub fn compute_tle_checksum(prefix_68: &str) -> u8 {
    compute_checksum(prefix_68)
}

pub(crate) fn compute_checksum(prefix_68: &str) -> u8 {
    let mut sum: u32 = 0;
    for b in prefix_68.bytes() {
        if b.is_ascii_digit() {
            sum += (b - b'0') as u32;
        } else if b == b'-' {
            sum += 1;
        }
    }
    (sum % 10) as u8
}

/// Parse a classic 2-line TLE.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::parse_tle;
/// let tle = parse_tle(
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// assert_eq!(tle.norad_id.0, 25544);
/// ```
pub fn parse_tle(line1: &str, line2: &str) -> Result<TLE, TleError> {
    parse_tle_inner(None, line1, line2)
}

/// Parse a 3-line TLE (name + classic 2-line).
///
/// The leading `"0 "` prefix is stripped from `name` if present (per the
/// 3LE convention used by Celestrak).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::parse_3le;
/// let tle = parse_3le(
///     "ISS (ZARYA)",
///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
/// ).unwrap();
/// assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
/// ```
pub fn parse_3le(name: &str, line1: &str, line2: &str) -> Result<TLE, TleError> {
    let trimmed = name.trim_start_matches("0 ").trim();
    let owned = if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_string())
    };
    parse_tle_inner(owned, line1, line2)
}

fn parse_tle_inner(name: Option<String>, line1: &str, line2: &str) -> Result<TLE, TleError> {
    if line1.len() != 69 {
        return Err(TleError::BadLength {
            line: 1,
            got: line1.len(),
        });
    }
    if line2.len() != 69 {
        return Err(TleError::BadLength {
            line: 2,
            got: line2.len(),
        });
    }
    let l1b = line1.as_bytes();
    let l2b = line2.as_bytes();
    if l1b[0] as char != '1' {
        return Err(TleError::BadLeadingChar {
            line: 1,
            expected: '1',
            found: l1b[0] as char,
        });
    }
    if l2b[0] as char != '2' {
        return Err(TleError::BadLeadingChar {
            line: 2,
            expected: '2',
            found: l2b[0] as char,
        });
    }
    validate_tle_checksum(line1)?;
    validate_tle_checksum(line2)?;

    let norad_l1 = SatelliteNumber::parse(&line1[2..7])?;
    let classification = Classification::from_char(line1.as_bytes()[7] as char)?;
    let intl_des = InternationalDesignator(line1[9..17].trim().to_string());
    let epoch_year_2 = parse_int_field(&line1[18..20], "epoch_year")?;
    let epoch_day = parse_float_field(&line1[20..32], "epoch_day")?;
    let mean_motion_dot = parse_signed_decimal(&line1[33..43], "mean_motion_dot")?;
    let mean_motion_ddot = parse_assumed_decimal_exponent(&line1[44..52], "mean_motion_ddot")?;
    let bstar = parse_assumed_decimal_exponent(&line1[53..61], "bstar")?;
    let element_set_number =
        parse_int_field(line1[64..68].trim_start(), "element_set_number")? as u16;

    let norad_l2 = SatelliteNumber::parse(&line2[2..7])?;
    if norad_l1 != norad_l2 {
        return Err(TleError::MismatchedSatelliteNumber {
            l1: norad_l1.0,
            l2: norad_l2.0,
        });
    }
    let inclination = Degrees::new(parse_float_field(&line2[8..16], "inclination")?);
    let raan = Degrees::new(parse_float_field(&line2[17..25], "raan")?);
    let eccentricity = {
        let raw = line2[26..33].trim();
        if raw.is_empty() || !raw.chars().all(|c| c.is_ascii_digit()) {
            return Err(TleError::InvalidNumber {
                field: "eccentricity",
                raw: raw.into(),
            });
        }
        let scaled: u64 = raw.parse().map_err(|_| TleError::InvalidNumber {
            field: "eccentricity",
            raw: raw.into(),
        })?;
        scaled as f64 * 10f64.powi(-(raw.len() as i32))
    };
    let argp = Degrees::new(parse_float_field(&line2[34..42], "argument_of_perigee")?);
    let mean_anomaly = Degrees::new(parse_float_field(&line2[43..51], "mean_anomaly")?);
    let mean_motion_revs_per_day = parse_float_field(&line2[52..63], "mean_motion")?;
    let mean_motion = AngularRate::<Turn, Day>::new(mean_motion_revs_per_day);
    let revolution_number_at_epoch =
        parse_int_field(line2[63..68].trim_start(), "revolution_number")? as u32;

    let epoch = epoch_from_year_doy(expand_two_digit_year(epoch_year_2 as i32), epoch_day)?;

    Ok(TLE {
        name,
        norad_id: norad_l1,
        classification,
        international_designator: intl_des,
        epoch,
        mean_motion_dot,
        mean_motion_ddot,
        bstar,
        element_set_number,
        revolution_number_at_epoch,
        inclination,
        raan,
        eccentricity,
        argument_of_perigee: argp,
        mean_anomaly,
        mean_motion,
    })
}

/// TLE 2-digit year convention: 57..99 → 1957..1999, 00..56 → 2000..2056.
pub(crate) fn expand_two_digit_year(yy: i32) -> i32 {
    if yy >= 57 {
        1900 + yy
    } else {
        2000 + yy
    }
}

pub(crate) fn epoch_from_year_doy(year: i32, day_of_year: f64) -> Result<Time<UTC>, TleError> {
    if !day_of_year.is_finite() || !(1.0..367.0).contains(&day_of_year) {
        return Err(TleError::InvalidEpoch {
            year,
            day_of_year,
            reason: "day-of-year out of range",
        });
    }
    let day_int = day_of_year.floor() as i64;
    let frac_day = day_of_year - day_int as f64;
    let date = NaiveDate::from_yo_opt(year, day_int as u32).ok_or(TleError::InvalidEpoch {
        year,
        day_of_year,
        reason: "calendar date does not exist",
    })?;
    let nanos = (frac_day * 86_400.0 * 1e9).round() as i64;
    let dt =
        Utc.from_utc_datetime(&date.and_hms_opt(0, 0, 0).unwrap()) + Duration::nanoseconds(nanos);
    Time::<UTC>::try_from_chrono(dt).map_err(|e| TleError::EpochConversion(format!("{e:?}")))
}

fn parse_int_field(raw: &str, field: &'static str) -> Result<i64, TleError> {
    let s = raw.trim();
    s.parse::<i64>().map_err(|_| TleError::InvalidNumber {
        field,
        raw: raw.into(),
    })
}

fn parse_float_field(raw: &str, field: &'static str) -> Result<f64, TleError> {
    let s = raw.trim();
    s.parse::<f64>().map_err(|_| TleError::InvalidNumber {
        field,
        raw: raw.into(),
    })
}

fn parse_signed_decimal(raw: &str, field: &'static str) -> Result<f64, TleError> {
    let s = raw.trim();
    s.parse::<f64>().map_err(|_| TleError::InvalidNumber {
        field,
        raw: raw.into(),
    })
}

/// Parse the TLE assumed-decimal-with-exponent encoding (e.g. `" 12345-3"` →
/// `0.12345e-3`).
pub(crate) fn parse_assumed_decimal_exponent(
    raw: &str,
    field: &'static str,
) -> Result<f64, TleError> {
    let s = raw.trim();
    if s.is_empty() {
        return Ok(0.0);
    }
    let (sign, rest) = match s.as_bytes()[0] {
        b'-' => (-1.0_f64, &s[1..]),
        b'+' => (1.0, &s[1..]),
        _ => (1.0, s),
    };
    if rest.len() < 2 {
        return Err(TleError::InvalidNumber {
            field,
            raw: raw.into(),
        });
    }
    let exp_sep = rest.rfind(['-', '+']).ok_or(TleError::InvalidNumber {
        field,
        raw: raw.into(),
    })?;
    let mantissa_digits = &rest[..exp_sep];
    let exp_str = &rest[exp_sep..];
    if !mantissa_digits.chars().all(|c| c.is_ascii_digit()) || mantissa_digits.is_empty() {
        return Err(TleError::InvalidNumber {
            field,
            raw: raw.into(),
        });
    }
    let mantissa: f64 = mantissa_digits
        .parse::<u64>()
        .map_err(|_| TleError::InvalidNumber {
            field,
            raw: raw.into(),
        })? as f64
        * 10f64.powi(-(mantissa_digits.len() as i32));
    let exp: i32 = exp_str.parse().map_err(|_| TleError::InvalidNumber {
        field,
        raw: raw.into(),
    })?;
    Ok(sign * mantissa * 10f64.powi(exp))
}
