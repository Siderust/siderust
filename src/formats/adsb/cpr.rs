// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Compact Position Reporting (CPR) decoder
//!
//! ## Scientific scope
//!
//! CPR is the position encoding used in ADS-B airborne position messages
//! (TC 9–18).  A single CPR-encoded position is ambiguous; geodetic
//! coordinates can be recovered by either:
//!
//! - **Global decoding** — combining one even and one odd frame transmitted
//!   within 10 seconds.
//! - **Local decoding** — using one frame together with a reference position
//!   (e.g. a known receiver location) accurate to within ≈180 NM.
//!
//! This implementation follows the algorithm in Sun (2021), §4.2 and §4.3,
//! which is equivalent to the algorithm in RTCA DO-260B §2.2.3.2.6.
//!
//! ## References
//!
//! - Sun, J. (2021). *The 1090 MHz Riddle*, 2nd ed. TU Delft Open. §4.2–4.3.
//! - RTCA (2009). *DO-260B*, §2.2.3.2.6.

use crate::formats::{FileLocation, FormatError};
use crate::qtty::Degrees;

/// Number of longitude zones at the equator.
const NZ: f64 = 15.0;

/// CPR latitude zone width for even (j=0) and odd (j=1) frames.
///
/// `dlat(j) = 360 / (4·NZ − j)`
#[inline]
fn dlat(odd: bool) -> f64 {
    360.0 / (4.0 * NZ - if odd { 1.0 } else { 0.0 })
}

/// Number of longitude zones at a given latitude.
///
/// Returns 1 when latitude is near the poles.
///
/// Based on DO-260B eq. 2-15 / Sun 2021 eq. 4.13.
fn nl(lat: f64) -> f64 {
    if lat.abs() >= 87.0 {
        return 1.0;
    }
    let cos_lat = lat.to_radians().cos();
    let numer = 1.0 - (std::f64::consts::PI / (2.0 * NZ)).cos();
    let arg = 1.0 - numer / (cos_lat * cos_lat);
    let n = (2.0 * std::f64::consts::PI / arg.acos()).floor();
    // Cap at 59: at lat=0 the formula evaluates to exactly 60.
    n.clamp(1.0, 59.0)
}

/// Decode a globally-unambiguous geodetic position from an even and odd
/// CPR frame pair.
///
/// Both frames must have been received within 10 seconds of each other.
/// The more recent frame determines which solution is selected.
///
/// # Arguments
///
/// - `lat0` / `lon0` — 17-bit CPR latitude / longitude from the **even** frame.
/// - `lat1` / `lon1` — 17-bit CPR latitude / longitude from the **odd** frame.
/// - `odd_is_latest` — `true` if the odd frame was received last.
///
/// # Returns
///
/// `Ok((Degrees, Degrees))` = (latitude, longitude) of the aircraft.
///
/// # Errors
///
/// [`FormatError`] when the even/odd pair is inconsistent (NL mismatch —
/// meaning the pair is too old and straddles a latitude zone boundary).
///
/// # Examples
///
/// ```rust
/// use siderust::formats::adsb::cpr::decode_globally;
///
/// // CPR values roundtrip-verified for ≈ 52.257°N, 3.919°E
/// let (lat, lon) = decode_globally(93000, 51372, 73974, 49949, false).unwrap();
/// assert!((lat.value() - 52.257).abs() < 0.01, "lat = {}", lat.value());
/// assert!((lon.value() - 3.919).abs() < 0.02, "lon = {}", lon.value());
/// ```
pub fn decode_globally(
    lat0: u32,
    lon0: u32,
    lat1: u32,
    lon1: u32,
    odd_is_latest: bool,
) -> Result<(Degrees, Degrees), FormatError> {
    let dlat0 = dlat(false);
    let dlat1 = dlat(true);

    let j = (59.0 * (lat0 as f64 / 131_072.0) - 60.0 * (lat1 as f64 / 131_072.0) + 0.5).floor();

    let lat_e = dlat0 * (j.rem_euclid(60.0) + lat0 as f64 / 131_072.0);
    let lat_o = dlat1 * (j.rem_euclid(59.0) + lat1 as f64 / 131_072.0);

    let lat_e = if lat_e >= 270.0 { lat_e - 360.0 } else { lat_e };
    let lat_o = if lat_o >= 270.0 { lat_o - 360.0 } else { lat_o };

    // NL check: both frames must be in the same latitude band
    if (nl(lat_e) - nl(lat_o)).abs() > 0.5 {
        return Err(FormatError::located(
            "ADS-B DO-260B §2.2.3.2.6",
            FileLocation::default(),
            "NL mismatch: even/odd frame pair straddles a latitude zone boundary",
        ));
    }

    let lat = if odd_is_latest { lat_o } else { lat_e };
    let nl_lat = nl(lat);

    let (lon, lon_ref) = if odd_is_latest {
        let dlon = if nl_lat <= 1.0 {
            360.0
        } else {
            360.0 / (nl_lat - 1.0)
        };
        let m = (lon0 as f64 * (nl_lat - 1.0) / 131_072.0 - lon1 as f64 * nl_lat / 131_072.0 + 0.5)
            .floor();
        (
            dlon * (m.rem_euclid(nl_lat - 1.0) + lon1 as f64 / 131_072.0),
            lat_o,
        )
    } else {
        let dlon = if nl_lat < 1.0 { 360.0 } else { 360.0 / nl_lat };
        let m = (lon0 as f64 * (nl_lat - 1.0) / 131_072.0 - lon1 as f64 * nl_lat / 131_072.0 + 0.5)
            .floor();
        (
            dlon * (m.rem_euclid(nl_lat) + lon0 as f64 / 131_072.0),
            lat_e,
        )
    };

    let _ = lon_ref; // used implicitly in lat selection above

    let lon = if lon >= 180.0 { lon - 360.0 } else { lon };

    Ok((Degrees::new(lat), Degrees::new(lon)))
}

/// Decode a locally-unambiguous geodetic position from a single CPR frame
/// and a reference position.
///
/// The reference position must be within ≈180 NM (≈333 km) of the aircraft.
///
/// # Arguments
///
/// - `cpr_lat` / `cpr_lon` — 17-bit CPR fields from the frame.
/// - `odd` — `true` if this is an odd frame.
/// - `ref_lat` / `ref_lon` — reference position in degrees.
///
/// # Returns
///
/// `Ok((Degrees, Degrees))` = (latitude, longitude).
///
/// # Examples
///
/// ```rust
/// use siderust::formats::adsb::cpr::decode_locally;
///
/// // Decode locally with reference near Amsterdam (52°N, 4.9°E)
/// let (lat, lon) = decode_locally(93000, 51372, false, 52.0, 4.9).unwrap();
/// assert!((lat.value() - 52.257).abs() < 0.01, "lat = {}", lat.value());
/// ```
pub fn decode_locally(
    cpr_lat: u32,
    cpr_lon: u32,
    odd: bool,
    ref_lat: f64,
    ref_lon: f64,
) -> Result<(Degrees, Degrees), FormatError> {
    let d_lat = dlat(odd);
    let j = (ref_lat / d_lat - cpr_lat as f64 / 131_072.0 + 0.5).floor();
    let lat = d_lat * (j + cpr_lat as f64 / 131_072.0);
    let lat = if lat >= 270.0 { lat - 360.0 } else { lat };

    let nl_val = nl(lat);
    let n_lon = if odd {
        (nl_val - 1.0).max(1.0)
    } else {
        nl_val.max(1.0)
    };
    let d_lon = 360.0 / n_lon;
    let m = (ref_lon / d_lon - cpr_lon as f64 / 131_072.0 + 0.5).floor();
    let lon = d_lon * (m + cpr_lon as f64 / 131_072.0);
    let lon = if lon >= 180.0 { lon - 360.0 } else { lon };

    Ok((Degrees::new(lat), Degrees::new(lon)))
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // CPR values roundtrip-verified for lat ≈ 52.257°N, lon ≈ 3.919°E
    // (extracted from Sun 2021 §4.2 frame 8D40621D58C382D690C8AC2863A7 + synthetic odd companion).
    // even: lat_cpr=93000 lon_cpr=51372; odd: lat_cpr=73974 lon_cpr=49949
    #[test]
    fn global_decode_sun_example() {
        let (lat, lon) = decode_globally(93_000, 51_372, 73_974, 49_949, false).unwrap();
        assert!((lat.value() - 52.257).abs() < 0.01, "lat={}", lat.value());
        assert!((lon.value() - 3.919).abs() < 0.02, "lon={}", lon.value());
    }

    #[test]
    fn local_decode_sun_example() {
        // Even frame, reference near Amsterdam (52°N, 4.9°E)
        let (lat, lon) = decode_locally(93_000, 51_372, false, 52.0, 4.9).unwrap();
        assert!((lat.value() - 52.257).abs() < 0.01, "lat={}", lat.value());
        assert!((lon.value() - 3.919).abs() < 0.1, "lon={}", lon.value());
    }

    #[test]
    fn nl_equator() {
        // At equator NL should be 59
        assert_eq!(nl(0.0) as u32, 59);
    }

    #[test]
    fn nl_polar() {
        // Above 87° NL = 1
        assert_eq!(nl(87.5) as u32, 1);
    }
}
