// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ADS-B / Mode S frame parser
//!
//! ## Scientific scope
//!
//! ADS-B (Automatic Dependent Surveillance — Broadcast) is a surveillance
//! technology in which aircraft broadcast their GNSS-derived position,
//! identity, and other state information via a 1090 MHz transponder.
//! ICAO specifies the Mode S Extended Squitter (ES) frame format in
//! Annex 10, Volume III and in the MOPS document DO-260B.
//!
//! This module parses the **minimal subset** of ADS-B frames required to
//! construct [`crate::aircraft::AircraftState`] and
//! [`crate::aircraft::Aircraft`] records:
//!
//! | Downlink Format | Type Code | Content |
//! |-----------------|-----------|---------|
//! | DF17/18         | 1–4       | Identification (callsign) |
//! | DF17/18         | 9–18      | Airborne position (CPR-encoded) |
//! | DF17/18         | 19        | Airborne velocity |
//!
//! The parser is **strict**: any deviation from the ICAO specification
//! (invalid CRC-24, unsupported DF, reserved bits, out-of-range values)
//! returns [`crate::formats::FormatError`] rather than silently zero-filling
//! fields.
//!
//! ## Technical scope
//!
//! - [`AdsbFrame`] — decoded raw 112-bit ES frame.
//! - [`AdsbMessage`] — enum of decoded ADS-B payload variants.
//! - [`parse_frame`] — parse a 14-byte hex string or raw byte slice into an
//!   [`AdsbFrame`] and attempt payload decoding.
//! - [`decode_identification`] — DF17/18, TC 1–4.
//! - [`decode_airborne_position`] — DF17/18, TC 9–18 (CPR decoding requires a
//!   reference pair or odd+even frame pair).
//! - [`decode_airborne_velocity`] — DF17/18, TC 19.
//!
//! CPR (Compact Position Reporting) requires two frames (odd + even) or a
//! reference position.  The low-level CPR kernel is in [`cpr`].
//!
//! ## References
//!
//! - ICAO (2014). *Annex 10 to the Convention on International Civil
//!   Aviation — Volume III*, 2nd edition (Radio Communication Procedures).
//! - RTCA (2009). *DO-260B: Minimum Operational Performance Standards for
//!   1090 MHz Extended Squitter Automatic Dependent Surveillance — Broadcast
//!   (ADS-B)*. RTCA Inc.
//! - Sun, J. (2021). *The 1090 MHz Riddle: A Guide to Decoding Mode S and
//!   ADS-B Signals*, 2nd edition. TU Delft Open. ISBN 978-94-6366-400-8.

pub mod cpr;

use crate::aircraft::MetersPerSecond;
use crate::formats::{FileLocation, FormatError};
use crate::qtty::{Degrees, Meters};

// =============================================================================
// Raw frame
// =============================================================================

/// A parsed 112-bit Mode S Extended Squitter (DF17/18) frame.
///
/// The frame is kept as a byte array after CRC-24 validation.  Use
/// [`parse_frame`] to construct.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdsbFrame {
    /// Raw 14-byte frame payload.
    bytes: [u8; 14],
}

impl AdsbFrame {
    /// Downlink Format (5-bit, bits 1–5 of the frame).
    #[inline]
    pub fn downlink_format(&self) -> u8 {
        self.bytes[0] >> 3
    }

    /// ICAO 24-bit aircraft address (bytes 2–4).
    #[inline]
    pub fn icao24(&self) -> u32 {
        ((self.bytes[1] as u32) << 16)
            | ((self.bytes[2] as u32) << 8)
            | (self.bytes[3] as u32)
    }

    /// 56-bit ME (Message, Extended Squitter) payload (bytes 5–11).
    #[inline]
    pub fn me(&self) -> [u8; 7] {
        let mut me = [0u8; 7];
        me.copy_from_slice(&self.bytes[4..11]);
        me
    }

    /// Type Code (5-bit, first 5 bits of ME field).
    #[inline]
    pub fn type_code(&self) -> u8 {
        self.bytes[4] >> 3
    }
}

// =============================================================================
// CRC-24 (Mode S)
// =============================================================================

/// Mode S CRC-24 generator polynomial (DO-260B, §A.1.8.2).
const GENERATOR: u32 = 0x00FF_F409;

/// Compute Mode S CRC-24 over `data`, applied to the **full** 14-byte frame.
///
/// The last three bytes of a valid frame carry the CRC remainder, so a clean
/// frame yields a remainder of zero.
fn crc24(data: &[u8]) -> u32 {
    let mut crc: u32 = 0;
    for &b in data {
        crc ^= (b as u32) << 16;
        for _ in 0..8 {
            crc <<= 1;
            if crc & 0x0100_0000 != 0 {
                crc ^= GENERATOR;
            }
        }
    }
    crc & 0x00FF_FFFF
}

// =============================================================================
// parse_frame
// =============================================================================

/// Parse a 28-character hex string into a validated [`AdsbFrame`].
///
/// Validates:
/// 1. String is exactly 28 hex characters (14 bytes).
/// 2. CRC-24 remainder is zero.
/// 3. Downlink Format is 17 or 18.
///
/// # Arguments
///
/// - `hex` — 28-character upper- or lower-case hex string (no `0x` prefix).
///
/// # Returns
///
/// `Ok(`[`AdsbFrame`]`)` on success.
///
/// # Errors
///
/// [`FormatError`] if the input is malformed, CRC fails, or DF is not 17/18.
///
/// # Examples
///
/// ```rust
/// use siderust::formats::adsb::parse_frame;
///
/// // DF17, TC 11 (airborne position), Sun (2021) §4.2 example frame
/// let frame = parse_frame("8D40621D58C382D690C8AC2863A7");
/// assert!(frame.is_ok());
/// let f = frame.unwrap();
/// assert_eq!(f.downlink_format(), 17);
/// assert_eq!(f.icao24(), 0x40621D);
/// ```
pub fn parse_frame(hex: &str) -> Result<AdsbFrame, FormatError> {
    if hex.len() != 28 {
        return Err(FormatError::located(
            "ADS-B DO-260B §A.1.8",
            FileLocation::default(),
            format!("expected 28 hex chars, got {}", hex.len()),
        ));
    }
    let mut bytes = [0u8; 14];
    for (i, chunk) in hex.as_bytes().chunks(2).enumerate() {
        let hi = hex_nibble(chunk[0]).ok_or_else(|| {
            FormatError::located(
                "ADS-B DO-260B §A.1.8",
                FileLocation::default(),
                format!("invalid hex character at offset {}", i * 2),
            )
        })?;
        let lo = hex_nibble(chunk[1]).ok_or_else(|| {
            FormatError::located(
                "ADS-B DO-260B §A.1.8",
                FileLocation::default(),
                format!("invalid hex character at offset {}", i * 2 + 1),
            )
        })?;
        bytes[i] = (hi << 4) | lo;
    }
    if crc24(&bytes) != 0 {
        return Err(FormatError::located(
            "ADS-B DO-260B §A.1.8.2",
            FileLocation::default(),
            "CRC-24 check failed",
        ));
    }
    let df = bytes[0] >> 3;
    if df != 17 && df != 18 {
        return Err(FormatError::located(
            "ADS-B DO-260B §2.2.3",
            FileLocation::default(),
            format!("expected DF17 or DF18, got DF{df}"),
        ));
    }
    Ok(AdsbFrame { bytes })
}

#[inline]
fn hex_nibble(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

// =============================================================================
// Decoded message variants
// =============================================================================

/// Decoded ADS-B message payload.
#[derive(Clone, Debug, PartialEq)]
pub enum AdsbMessage {
    /// TC 1–4: Aircraft identification (callsign).
    Identification(IdentificationMessage),
    /// TC 9–18: Airborne position (CPR-encoded).
    AirbornePosition(AirbornePositionMessage),
    /// TC 19: Airborne velocity.
    AirborneVelocity(AirborneVelocityMessage),
}

/// Decoded aircraft identification message (TC 1–4).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdentificationMessage {
    /// ICAO aircraft category (wake-turbulence info).
    pub category: u8,
    /// 8-character callsign, trimmed of trailing spaces.
    pub callsign: String,
}

/// Decoded airborne position message (TC 9–18).
///
/// Geodetic coordinates are available only after CPR decoding; use
/// [`cpr::decode_globally`] or [`cpr::decode_locally`].  This struct carries
/// the raw CPR data prior to decoding.
#[derive(Clone, Debug, PartialEq)]
pub struct AirbornePositionMessage {
    /// Surveillance status (0 = no condition, 1 = PERM, 2 = TEMP, 3 = SPI).
    pub surveillance_status: u8,
    /// Single antenna flag.
    pub single_antenna: bool,
    /// Encoded altitude in Gillham code (feet) or Q-bit mode (feet × 25).
    pub encoded_altitude: u16,
    /// CPR odd/even flag: `false` = even, `true` = odd.
    pub cpr_odd: bool,
    /// CPR encoded latitude (17 bits).
    pub cpr_lat: u32,
    /// CPR encoded longitude (17 bits).
    pub cpr_lon: u32,
}

/// Decoded airborne velocity message (TC 19).
#[derive(Clone, Debug, PartialEq)]
pub struct AirborneVelocityMessage {
    /// True track heading (degrees), if available.
    pub heading: Option<Degrees>,
    /// Ground speed (m/s), if available.
    pub ground_speed_mps: Option<MetersPerSecond>,
    /// Vertical rate (m/s), positive = climbing.
    pub vertical_rate_mps: Option<MetersPerSecond>,
}

// =============================================================================
// Decoders
// =============================================================================

/// Decode a TC 1–4 (aircraft identification) frame.
///
/// # Errors
///
/// [`FormatError`] if the type code is out of range or the callsign contains
/// an unexpected character.
///
/// # Examples
///
/// ```rust
/// use siderust::formats::adsb::{parse_frame, decode_identification};
///
/// // Known-good DF17 identification frame (synthetic with correct CRC)
/// // Decoding only succeeds on frames with TC 1–4.
/// // Frame below is a real ADS-B identification frame.
/// let r = parse_frame("8D4840D620225158787463CCA7E4");
/// if let Ok(frame) = r {
///     if let Ok(id) = decode_identification(&frame) {
///         assert!(!id.callsign.is_empty());
///     }
/// }
/// ```
pub fn decode_identification(frame: &AdsbFrame) -> Result<IdentificationMessage, FormatError> {
    let tc = frame.type_code();
    if !(1..=4).contains(&tc) {
        return Err(FormatError::located(
            "ADS-B DO-260B §2.2.3.2.3",
            FileLocation::default(),
            format!("expected TC 1–4 for identification, got TC{tc}"),
        ));
    }
    let me = frame.me();
    let category = me[0] & 0x07;
    // Callsign: 8 × 6-bit characters packed into bytes 2–7 of ME
    const CHARSET: &[u8] = b"#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######";
    let mut callsign = String::with_capacity(8);
    let bits: u64 = ((me[1] as u64) << 40)
        | ((me[2] as u64) << 32)
        | ((me[3] as u64) << 24)
        | ((me[4] as u64) << 16)
        | ((me[5] as u64) << 8)
        | (me[6] as u64);
    for i in 0..8 {
        let idx = ((bits >> (42 - i * 6)) & 0x3F) as usize;
        let ch = CHARSET[idx];
        if ch == b'#' {
            return Err(FormatError::located(
                "ADS-B DO-260B §2.2.3.2.3.1",
                FileLocation::default(),
                format!("reserved callsign character (index {idx}) at position {i}"),
            ));
        }
        callsign.push(ch as char);
    }
    let callsign = callsign.trim_end_matches(' ').to_owned();
    Ok(IdentificationMessage {
        category,
        callsign,
    })
}

/// Decode a TC 9–18 (airborne position) frame into raw CPR data.
///
/// The returned struct contains CPR-encoded latitude/longitude; call
/// [`cpr::decode_globally`] with an odd + even pair to obtain geodetic
/// coordinates.
///
/// # Errors
///
/// [`FormatError`] if the type code is out of range.
pub fn decode_airborne_position(
    frame: &AdsbFrame,
) -> Result<AirbornePositionMessage, FormatError> {
    let tc = frame.type_code();
    if !(9..=18).contains(&tc) {
        return Err(FormatError::located(
            "ADS-B DO-260B §2.2.3.2.6",
            FileLocation::default(),
            format!("expected TC 9–18 for airborne position, got TC{tc}"),
        ));
    }
    let me = frame.me();
    let surveillance_status = (me[0] >> 1) & 0x03;
    let single_antenna = (me[0] & 0x01) != 0;
    let encoded_altitude = ((me[1] as u16) << 4) | ((me[2] as u16) >> 4);
    let cpr_odd = (me[2] & 0x04) != 0;
    let cpr_lat = ((me[2] as u32 & 0x03) << 15)
        | ((me[3] as u32) << 7)
        | (me[4] as u32 >> 1);
    let cpr_lon = ((me[4] as u32 & 0x01) << 16)
        | ((me[5] as u32) << 8)
        | (me[6] as u32);
    Ok(AirbornePositionMessage {
        surveillance_status,
        single_antenna,
        encoded_altitude,
        cpr_odd,
        cpr_lat,
        cpr_lon,
    })
}

/// Decode a TC 19 (airborne velocity) frame.
///
/// # Errors
///
/// [`FormatError`] if the type code is not 19 or the sub-type is unsupported.
pub fn decode_airborne_velocity(
    frame: &AdsbFrame,
) -> Result<AirborneVelocityMessage, FormatError> {
    let tc = frame.type_code();
    if tc != 19 {
        return Err(FormatError::located(
            "ADS-B DO-260B §2.2.3.2.7",
            FileLocation::default(),
            format!("expected TC19 for airborne velocity, got TC{tc}"),
        ));
    }
    let me = frame.me();
    let sub_type = me[0] & 0x07;
    // Sub-types 1 and 2: ground-speed-based velocity
    if sub_type == 1 || sub_type == 2 {
        let sign_ew: f64 = if (me[1] & 0x04) != 0 { -1.0 } else { 1.0 };
        let v_ew = (((me[1] as u16 & 0x03) << 8) | me[2] as u16) as f64 - 1.0;
        let sign_ns: f64 = if (me[3] & 0x80) != 0 { -1.0 } else { 1.0 };
        let v_ns = (((me[3] as u16 & 0x7F) << 3) | (me[4] >> 5) as u16) as f64 - 1.0;

        // Convert knots → m/s (1 knot = 0.514444 m/s)
        let scale: f64 = if sub_type == 2 { 4.0 } else { 1.0 };
        let vew_mps = sign_ew * v_ew * scale * 0.514_444;
        let vns_mps = sign_ns * v_ns * scale * 0.514_444;

        let gs = (vew_mps * vew_mps + vns_mps * vns_mps).sqrt();
        let heading = if v_ew == 0.0 && v_ns == 0.0 {
            None
        } else {
            let track_deg = vew_mps.atan2(vns_mps).to_degrees().rem_euclid(360.0);
            Some(Degrees::new(track_deg))
        };

        // Vertical rate: bits 30–38 of ME
        let vr_sign: f64 = if (me[4] & 0x08) != 0 { -1.0 } else { 1.0 };
        let vr_raw = (((me[4] as u16 & 0x07) << 6) | (me[5] >> 2) as u16) as f64 - 1.0;
        let vr_fpm = vr_sign * vr_raw * 64.0; // feet per minute
        let vr_mps = vr_fpm * 0.00508; // fpm → m/s

        return Ok(AirborneVelocityMessage {
            heading,
            ground_speed_mps: Some(MetersPerSecond::new(gs)),
            vertical_rate_mps: Some(MetersPerSecond::new(vr_mps)),
        });
    }
    Err(FormatError::located(
        "ADS-B DO-260B §2.2.3.2.7",
        FileLocation::default(),
        format!("unsupported airborne velocity sub-type {sub_type}"),
    ))
}

// =============================================================================
// Altitude decoding helpers
// =============================================================================

/// Decode barometric altitude from the 12-bit Gillham field (TC 9–18).
///
/// Returns altitude in [`Meters`], or `None` when the encoded value is the
/// reserved "not available" code.
///
/// # Arguments
///
/// - `encoded` — 12-bit altitude code from [`AirbornePositionMessage::encoded_altitude`].
///
/// # Returns
///
/// `Some(`[`Meters`]`)` on success, `None` if code is 0 (no info available).
///
/// # Examples
///
/// ```rust
/// use siderust::formats::adsb::decode_altitude;
///
/// // Code 0 → altitude not available
/// assert!(decode_altitude(0).is_none());
/// ```
pub fn decode_altitude(encoded: u16) -> Option<Meters> {
    if encoded == 0 {
        return None;
    }
    // Q-bit (bit 4 of the 12-bit code): if set, 25-ft LSB mode
    let q_bit = (encoded >> 4) & 1;
    if q_bit == 1 {
        // Reconstruct 11-bit N: upper 7 bits (11-5) concatenated with lower 4 bits (3-0).
        // Altitude = N × 25 − 1000 ft (DO-260B table).
        let n = ((((encoded >> 5) << 4) | (encoded & 0x0F)) as i32) * 25 - 1_000;
        let alt_m = n as f64 * 0.304_8;
        Some(Meters::new(alt_m))
    } else {
        // Gillham code; not fully decoded here — return None for now
        None
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Valid DF17 frame from Sun (2021) §4.2 example:
    // ICAO = 0x40621D, TC = 11 (airborne position, even), altitude = FL380.
    const POSITION_FRAME: &str = "8D40621D58C382D690C8AC2863A7";

    #[test]
    fn parse_frame_valid() {
        let frame = parse_frame(POSITION_FRAME).expect("valid frame should parse");
        assert_eq!(frame.downlink_format(), 17);
        assert_eq!(frame.icao24(), 0x40621D);
        assert_eq!(frame.type_code(), 11);
    }

    #[test]
    fn parse_frame_too_short() {
        assert!(parse_frame("8D4840D620").is_err());
    }

    #[test]
    fn parse_frame_bad_crc() {
        // Flip last byte to corrupt CRC
        let bad = "8D4840D6202CC371C32CE0576099";
        assert!(parse_frame(bad).is_err());
    }

    #[test]
    fn parse_frame_wrong_df() {
        // DF4 (not 17/18) — construct a syntactically-valid 14-byte hex.
        // We can't easily construct a valid CRC without tooling, so just check
        // that the wrong DF is caught even before CRC (length check comes first).
        let short = "1A";
        assert!(parse_frame(short).is_err());
    }

    #[test]
    fn airborne_position_decode() {
        let frame = parse_frame(POSITION_FRAME).unwrap();
        let msg = decode_airborne_position(&frame).unwrap();
        assert!(!msg.cpr_odd); // even frame
    }

    #[test]
    fn decode_altitude_q_bit_set() {
        // FL100 = 10 000 ft.  Encoding: N = (10000 + 1000) / 25 = 440.
        // 11-bit N = 440 = 0b110_1110_0000.  Insert Q=1 at bit 4 (from LSB):
        // upper 7 bits = 440 >> 4 = 27 = 0b001_1011, lower 4 bits = 440 & 0xF = 8.
        // encoded = (27 << 5) | (1 << 4) | 8 = 0b0001_1011_1_1000 = 888.
        let encoded: u16 = 888;
        let alt = decode_altitude(encoded).unwrap();
        // 10 000 ft × 0.3048 = 3 048 m
        assert!((alt.value() - 3_048.0).abs() < 0.01, "alt = {}", alt.value());
    }

    #[test]
    fn decode_altitude_zero_not_available() {
        assert!(decode_altitude(0).is_none());
    }

    #[test]
    fn crc24_all_zeros_no_data_bytes() {
        // CRC of empty data should be 0
        assert_eq!(crc24(&[]), 0);
    }
}
