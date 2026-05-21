//! # ANTEX antenna phase-centre offsets
//!
//! ## Scientific scope
//!
//! ANTEX files describe antenna phase-centre offsets and variations for
//! GNSS tracking equipment. This module focuses on the offset portion
//! of that standard because the north-east-up offsets are what the
//! GNSS observation pipeline needs to relate receiver/satellite mass
//! centres to the antenna phase centre.
//!
//! The current regime is intentionally narrow: phase-centre variation grids
//! are ignored, so this parser is suitable for offset-driven geometry
//! corrections but not for full antenna pattern modelling.
//!
//! ## Technical scope
//!
//! The main entry point is `read_antex`, which returns a nested catalog
//! keyed by antenna and frequency identifiers. Offsets are exposed through
//! `Pco` displacement vectors in the millimetre convention carried by the
//! file. ANTEX orders the components as north, east, up; this module maps that
//! ordering to `x`, `y`, `z`.
//!
//! This module only parses the exchange format; it does not decide how
//! those offsets are applied inside any specific observation model.
//!
//! ## References
//!
//! - IGS MGEX/IGS. (2015). ANTEX: The Antenna Exchange Format, Version 1.4.
//! - Montenbruck, O., Steigenberger, P., & Hauschild, A. (2015). Broadcast
//!   versus precise ephemerides for GNSS orbit determination. GPS
//!   Solutions, 19(2), 321-330.
use super::FormatError;
use affn::cartesian::Displacement;
use affn::frames::ReferenceFrame;
use qtty::length::Millimeter;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

/// ANTEX local north-east-up component frame.
///
/// Component mapping for Cartesian vectors in this frame is:
///
/// - `x`: north
/// - `y`: east
/// - `z`: up
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AntexNeu;

impl ReferenceFrame for AntexNeu {
    fn frame_name() -> &'static str {
        "ANTEX NEU"
    }
}

/// PCO entry for one antenna and one frequency.
///
/// ANTEX stores phase-centre offsets as north/east/up millimetres. This type
/// is a free displacement vector with `x/y/z = north/east/up`.
pub type Pco = Displacement<AntexNeu, Millimeter>;

/// Per-antenna PCO data, keyed by frequency identifier (e.g. "G01" for GPS L1).
pub type AntennaPco = HashMap<String, Pco>;

/// Parsed ANTEX file: `antenna_name -> per-frequency PCO`.
pub type AntexCatalog = HashMap<String, AntennaPco>;

/// Parse an ANTEX file (subset).
pub fn read_antex<R: Read>(rdr: R) -> Result<AntexCatalog, FormatError> {
    let br = BufReader::new(rdr);
    let mut catalog = AntexCatalog::new();
    let mut current_antenna: Option<String> = None;
    let mut current_pcos: AntennaPco = AntennaPco::new();
    let mut current_freq: Option<String> = None;

    for line in br.lines() {
        let line = line?;
        let label = line.get(60..).unwrap_or("").trim();
        let body = line.get(..60).unwrap_or("");
        match label {
            "START OF ANTENNA" => {
                current_antenna = None;
                current_pcos.clear();
                current_freq = None;
            }
            "TYPE / SERIAL NO" => {
                current_antenna = Some(body.trim().to_string());
            }
            "START OF FREQUENCY" => {
                current_freq = body.split_whitespace().next().map(|s| s.to_string());
            }
            "NORTH / EAST / UP" => {
                if let Some(freq) = &current_freq {
                    let parts: Vec<f64> = body
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();
                    if parts.len() == 3 {
                        current_pcos.insert(freq.clone(), Pco::new(parts[0], parts[1], parts[2]));
                    }
                }
            }
            "END OF FREQUENCY" => {
                current_freq = None;
            }
            "END OF ANTENNA" => {
                if let Some(name) = current_antenna.take() {
                    catalog.insert(name, std::mem::take(&mut current_pcos));
                }
            }
            _ => {}
        }
    }
    Ok(catalog)
}

fn header_line<W: Write>(w: &mut W, body: &str, label: &str) -> std::io::Result<()> {
    let body = if body.len() > 60 { &body[..60] } else { body };
    writeln!(w, "{:<60}{}", body, label)
}

/// Write an ANTEX 1.4 catalog matching the layout consumed by [`read_antex`].
///
/// Output is canonicalised: antennas are emitted in lexicographic order by
/// type/serial string, and frequencies within each antenna are sorted by
/// frequency identifier. Phase-centre **variation grids** are not written
/// (only the offset block is).
///
/// Round-trip: `read_antex(write_antex(cat))` reproduces `cat` exactly
/// when only PCO offsets are present.
///
/// # Errors
///
/// Returns [`FormatError::Io`] on write failures.
///
/// # Examples
///
/// ```
/// use siderust::formats::igs::antex::{read_antex, write_antex};
///
/// let src = "\
///                                                             START OF ANTENNA\n\
/// TINY-ANT  001                                               TYPE / SERIAL NO\n\
///    G01                                                      START OF FREQUENCY\n\
///       1.00      0.00     90.00                              NORTH / EAST / UP\n\
///                                                             END OF FREQUENCY\n\
///                                                             END OF ANTENNA\n";
/// let cat = read_antex(src.as_bytes()).unwrap();
/// let mut buf = Vec::new();
/// write_antex(&mut buf, &cat).unwrap();
/// let cat2 = read_antex(&buf[..]).unwrap();
/// assert_eq!(cat.len(), cat2.len());
/// ```
pub fn write_antex<W: Write>(w: &mut W, catalog: &AntexCatalog) -> Result<(), FormatError> {
    header_line(
        w,
        "     1.4            G                  ",
        "ANTEX VERSION / SYST",
    )?;
    header_line(
        w,
        "ANTEX                                                       ",
        "PCV TYPE / REFANT",
    )?;
    header_line(w, "", "END OF HEADER")?;

    let mut antennas: Vec<&String> = catalog.keys().collect();
    antennas.sort();
    for ant in antennas {
        header_line(w, "", "START OF ANTENNA")?;
        header_line(w, ant, "TYPE / SERIAL NO")?;
        let pcos = &catalog[ant];
        let mut freqs: Vec<&String> = pcos.keys().collect();
        freqs.sort();
        for freq in freqs {
            let pco = &pcos[freq];
            header_line(w, &format!("   {:<3}", freq), "START OF FREQUENCY")?;
            let body = format!(
                "{:10.2}{:10.2}{:10.2}",
                pco.x().value(),
                pco.y().value(),
                pco.z().value()
            );
            header_line(w, &body, "NORTH / EAST / UP")?;
            header_line(w, &format!("   {:<3}", freq), "END OF FREQUENCY")?;
        }
        header_line(w, "", "END OF ANTENNA")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
                                                            START OF ANTENNA
TEST-ANTENNA   ABC                                          TYPE / SERIAL NO
   G01                                                      START OF FREQUENCY
      1.50      0.20     90.00                              NORTH / EAST / UP
                                                            END OF FREQUENCY
   G02                                                      START OF FREQUENCY
      0.50     -0.30     85.00                              NORTH / EAST / UP
                                                            END OF FREQUENCY
                                                            END OF ANTENNA
";

    #[test]
    fn parses_pco_blocks() {
        let cat = read_antex(SAMPLE.as_bytes()).unwrap();
        let ant = &cat["TEST-ANTENNA   ABC"];
        let g01 = ant["G01"];
        assert!((g01.z().value() - 90.0).abs() < 1e-9);
        let g02 = ant["G02"];
        assert!((g02.y().value() + 0.3).abs() < 1e-9);
    }

    #[test]
    fn round_trips_through_writer() {
        let cat = read_antex(SAMPLE.as_bytes()).unwrap();
        let mut buf = Vec::new();
        write_antex(&mut buf, &cat).unwrap();
        let cat2 = read_antex(&buf[..]).unwrap();
        assert_eq!(cat.len(), cat2.len());
        for (name, pcos) in &cat {
            let pcos2 = cat2.get(name).expect("antenna missing after round-trip");
            assert_eq!(pcos.len(), pcos2.len());
            for (freq, pco) in pcos {
                let pco2 = pcos2[freq];
                assert!((pco.x().value() - pco2.x().value()).abs() < 1e-3);
                assert!((pco.y().value() - pco2.y().value()).abs() < 1e-3);
                assert!((pco.z().value() - pco2.z().value()).abs() < 1e-3);
            }
        }
    }
}
