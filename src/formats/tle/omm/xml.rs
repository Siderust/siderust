// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Minimal XML encoding for OMM (CCSDS 502.0 Annex E).
//!
//! This is a hand-rolled, schema-restricted parser/writer that handles the
//! flat OMM element set required for SGP4 round-trip. It deliberately does
//! not pull in a general XML library — see crate-level rationale in
//! `lib.rs`.

use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;

use super::{format_epoch, parse_epoch, Omm};
use crate::formats::tle::TleError;
use crate::formats::tle::{Classification, SatelliteNumber};

/// Render an [`Omm`] to its XML form (`<omm>...</omm>`).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{omm::{Omm, xml}, parse_3le};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let omm = Omm::from_tle(&parse_3le("ISS (ZARYA)", l1, l2).unwrap());
/// let xml = xml::write(&omm).unwrap();
/// assert!(xml.contains("<NORAD_CAT_ID>25544</NORAD_CAT_ID>"));
/// ```
pub fn write(omm: &Omm) -> Result<String, TleError> {
    let epoch = format_epoch(omm.epoch)?;
    Ok(format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <omm id=\"CCSDS_OMM_VERS\" version=\"2.0\">\n\
         <body><segment><metadata>\n\
         <OBJECT_NAME>{name}</OBJECT_NAME>\n\
         <OBJECT_ID>{id}</OBJECT_ID>\n\
         <CENTER_NAME>EARTH</CENTER_NAME>\n\
         <REF_FRAME>TEME</REF_FRAME>\n\
         <TIME_SYSTEM>UTC</TIME_SYSTEM>\n\
         <MEAN_ELEMENT_THEORY>SGP4</MEAN_ELEMENT_THEORY>\n\
         </metadata><data>\n\
         <meanElements>\n\
         <EPOCH>{epoch}</EPOCH>\n\
         <MEAN_MOTION>{mm:.8}</MEAN_MOTION>\n\
         <ECCENTRICITY>{ecc:.7}</ECCENTRICITY>\n\
         <INCLINATION>{inc:.4}</INCLINATION>\n\
         <RA_OF_ASC_NODE>{raan:.4}</RA_OF_ASC_NODE>\n\
         <ARG_OF_PERICENTER>{argp:.4}</ARG_OF_PERICENTER>\n\
         <MEAN_ANOMALY>{ma:.4}</MEAN_ANOMALY>\n\
         </meanElements>\n\
         <tleParameters>\n\
         <EPHEMERIS_TYPE>{eph}</EPHEMERIS_TYPE>\n\
         <CLASSIFICATION_TYPE>{cls}</CLASSIFICATION_TYPE>\n\
         <NORAD_CAT_ID>{cat}</NORAD_CAT_ID>\n\
         <ELEMENT_SET_NO>{esn}</ELEMENT_SET_NO>\n\
         <REV_AT_EPOCH>{rev}</REV_AT_EPOCH>\n\
         <BSTAR>{bstar:.8e}</BSTAR>\n\
         <MEAN_MOTION_DOT>{ndot:.8e}</MEAN_MOTION_DOT>\n\
         <MEAN_MOTION_DDOT>{nddot:.8e}</MEAN_MOTION_DDOT>\n\
         </tleParameters>\n\
         </data></segment></body></omm>\n",
        name = xml_escape(&omm.object_name),
        id = xml_escape(&omm.object_id),
        mm = omm.mean_motion.value(),
        ecc = omm.eccentricity,
        inc = omm.inclination.value(),
        raan = omm.ra_of_asc_node.value(),
        argp = omm.arg_of_pericenter.value(),
        ma = omm.mean_anomaly.value(),
        eph = omm.ephemeris_type,
        cls = omm.classification.as_char(),
        cat = omm.norad_id.0,
        esn = omm.element_set_no,
        rev = omm.rev_at_epoch,
        bstar = omm.bstar,
        ndot = omm.mean_motion_dot,
        nddot = omm.mean_motion_ddot,
    ))
}

/// Parse an OMM XML document.
///
/// Tolerant of whitespace, attribute order, and additional unknown
/// elements; strict on field presence and numeric format.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::omm::{xml, Omm, kvn};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let omm = Omm::from_tle(&siderust::formats::tle::parse_3le("ISS (ZARYA)", l1, l2).unwrap());
/// let s = xml::write(&omm).unwrap();
/// let omm2 = xml::read(&s).unwrap();
/// assert_eq!(omm2.norad_id, omm.norad_id);
/// ```
pub fn read(input: &str) -> Result<Omm, TleError> {
    use std::collections::HashMap;
    let mut kv: HashMap<String, String> = HashMap::new();
    for tag in TLE_TAGS {
        if let Some(v) = extract_tag(input, tag) {
            kv.insert((*tag).to_string(), v);
        }
    }

    let theory = kv
        .get("MEAN_ELEMENT_THEORY")
        .map(String::as_str)
        .unwrap_or("SGP4");
    if !theory.eq_ignore_ascii_case("SGP4") {
        return Err(TleError::OmmUnsupportedTheory(theory.to_string()));
    }

    let get = |k: &'static str| -> Result<&String, TleError> {
        kv.get(k).ok_or(TleError::OmmMissingField(k))
    };
    let getf = |k: &'static str| -> Result<f64, TleError> {
        let v = get(k)?;
        v.parse::<f64>().map_err(|_| TleError::OmmMalformed {
            format: "XML",
            reason: format!("<{k}> = {v:?} is not a float"),
        })
    };

    let object_name = xml_unescape(get("OBJECT_NAME")?);
    let object_id = xml_unescape(get("OBJECT_ID")?);
    let epoch = parse_epoch(get("EPOCH")?)?;
    let cls_raw = kv
        .get("CLASSIFICATION_TYPE")
        .map(String::as_str)
        .unwrap_or("U");
    let cls_char = cls_raw.chars().next().unwrap_or('U');
    let classification = Classification::from_char(cls_char)?;
    let norad_id = SatelliteNumber::parse(get("NORAD_CAT_ID")?)?;

    Ok(Omm {
        object_name,
        object_id,
        epoch,
        mean_motion: AngularRate::<Turn, Day>::new(getf("MEAN_MOTION")?),
        eccentricity: getf("ECCENTRICITY")?,
        inclination: Degrees::new(getf("INCLINATION")?),
        ra_of_asc_node: Degrees::new(getf("RA_OF_ASC_NODE")?),
        arg_of_pericenter: Degrees::new(getf("ARG_OF_PERICENTER")?),
        mean_anomaly: Degrees::new(getf("MEAN_ANOMALY")?),
        ephemeris_type: kv
            .get("EPHEMERIS_TYPE")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        classification,
        norad_id,
        element_set_no: kv
            .get("ELEMENT_SET_NO")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        rev_at_epoch: kv
            .get("REV_AT_EPOCH")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        bstar: getf("BSTAR").unwrap_or(0.0),
        mean_motion_dot: getf("MEAN_MOTION_DOT").unwrap_or(0.0),
        mean_motion_ddot: getf("MEAN_MOTION_DDOT").unwrap_or(0.0),
    })
}

const TLE_TAGS: &[&str] = &[
    "OBJECT_NAME",
    "OBJECT_ID",
    "CENTER_NAME",
    "REF_FRAME",
    "TIME_SYSTEM",
    "MEAN_ELEMENT_THEORY",
    "EPOCH",
    "MEAN_MOTION",
    "ECCENTRICITY",
    "INCLINATION",
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY",
    "EPHEMERIS_TYPE",
    "CLASSIFICATION_TYPE",
    "NORAD_CAT_ID",
    "ELEMENT_SET_NO",
    "REV_AT_EPOCH",
    "BSTAR",
    "MEAN_MOTION_DOT",
    "MEAN_MOTION_DDOT",
];

fn extract_tag(input: &str, tag: &str) -> Option<String> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let i = input.find(&open)?;
    let after = &input[i + open.len()..];
    let j = after.find(&close)?;
    Some(after[..j].trim().to_string())
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn xml_unescape(s: &str) -> String {
    s.replace("&apos;", "'")
        .replace("&quot;", "\"")
        .replace("&gt;", ">")
        .replace("&lt;", "<")
        .replace("&amp;", "&")
}
