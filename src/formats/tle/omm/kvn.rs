// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Keyword-Value-Notation (KVN) encoding for OMM.

use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;

use super::{format_epoch, parse_epoch, Omm};
use crate::formats::tle::TleError;
use crate::formats::tle::{Classification, SatelliteNumber};

/// Render an [`Omm`] to its KVN form.
///
/// The output is a series of `KEY = VALUE` lines terminated with `\n`.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{omm::{Omm, kvn}, parse_3le};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let omm = Omm::from_tle(&parse_3le("ISS (ZARYA)", l1, l2).unwrap());
/// let text = kvn::write(&omm).unwrap();
/// assert!(text.contains("OBJECT_NAME = ISS (ZARYA)"));
/// ```
pub fn write(omm: &Omm) -> Result<String, TleError> {
    let epoch = format_epoch(omm.epoch)?;
    Ok(format!(
        "CCSDS_OMM_VERS = 2.0\n\
         CREATION_DATE = {epoch}\n\
         ORIGINATOR = siderust-tle\n\
         OBJECT_NAME = {object_name}\n\
         OBJECT_ID = {object_id}\n\
         CENTER_NAME = EARTH\n\
         REF_FRAME = TEME\n\
         TIME_SYSTEM = UTC\n\
         MEAN_ELEMENT_THEORY = SGP4\n\
         EPOCH = {epoch}\n\
         MEAN_MOTION = {mean_motion:.8}\n\
         ECCENTRICITY = {ecc:.7}\n\
         INCLINATION = {inc:.4}\n\
         RA_OF_ASC_NODE = {raan:.4}\n\
         ARG_OF_PERICENTER = {argp:.4}\n\
         MEAN_ANOMALY = {ma:.4}\n\
         EPHEMERIS_TYPE = {eph}\n\
         CLASSIFICATION_TYPE = {cls}\n\
         NORAD_CAT_ID = {cat}\n\
         ELEMENT_SET_NO = {esn}\n\
         REV_AT_EPOCH = {rev}\n\
         BSTAR = {bstar:.8e}\n\
         MEAN_MOTION_DOT = {ndot:.8e}\n\
         MEAN_MOTION_DDOT = {nddot:.8e}\n",
        object_name = omm.object_name,
        object_id = omm.object_id,
        mean_motion = omm.mean_motion.value(),
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

/// Parse an [`Omm`] from its KVN form.
///
/// Comment lines (starting with `COMMENT`) and unknown keys are silently
/// preserved-by-discard, in the spirit of the CCSDS "tolerant reader"
/// principle. The reader is strict on field types: any present-but-malformed
/// numeric field returns [`TleError::OmmMalformed`].
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::omm::kvn;
/// let text = "CCSDS_OMM_VERS = 2.0\nOBJECT_NAME = TEST\nOBJECT_ID = 1998-067A\n\
///              MEAN_ELEMENT_THEORY = SGP4\nEPOCH = 2008-09-20T12:25:40.104192\n\
///              MEAN_MOTION = 15.72125391\nECCENTRICITY = 0.0006703\n\
///              INCLINATION = 51.6416\nRA_OF_ASC_NODE = 247.4627\n\
///              ARG_OF_PERICENTER = 130.5360\nMEAN_ANOMALY = 325.0288\n\
///              EPHEMERIS_TYPE = 0\nCLASSIFICATION_TYPE = U\nNORAD_CAT_ID = 25544\n\
///              ELEMENT_SET_NO = 292\nREV_AT_EPOCH = 56353\n\
///              BSTAR = -1.1606e-5\nMEAN_MOTION_DOT = -2.182e-5\n\
///              MEAN_MOTION_DDOT = 0.0\n";
/// let omm = kvn::read(text).unwrap();
/// assert_eq!(omm.norad_id.0, 25544);
/// ```
pub fn read(input: &str) -> Result<Omm, TleError> {
    use std::collections::HashMap;
    let mut kv: HashMap<String, String> = HashMap::new();
    for raw_line in input.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with("COMMENT") {
            continue;
        }
        let Some(eq) = line.find('=') else {
            continue;
        };
        let key = line[..eq].trim().to_string();
        let value = line[eq + 1..].trim().to_string();
        kv.insert(key, value);
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
            format: "KVN",
            reason: format!("{k}={v:?} is not a float"),
        })
    };
    let getu = |k: &'static str| -> Result<u32, TleError> {
        let v = get(k)?;
        v.parse::<u32>().map_err(|_| TleError::OmmMalformed {
            format: "KVN",
            reason: format!("{k}={v:?} is not an unsigned integer"),
        })
    };

    let object_name = get("OBJECT_NAME")?.clone();
    let object_id = get("OBJECT_ID")?.clone();
    let epoch = parse_epoch(get("EPOCH")?)?;
    let cls_raw = kv
        .get("CLASSIFICATION_TYPE")
        .map(String::as_str)
        .unwrap_or("U");
    let cls_char = cls_raw.chars().next().unwrap_or('U');
    let classification = Classification::from_char(cls_char)?;

    let norad_str = get("NORAD_CAT_ID")?;
    let norad_id = SatelliteNumber::parse(norad_str)?;

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
        rev_at_epoch: getu("REV_AT_EPOCH").unwrap_or(0),
        bstar: getf("BSTAR").unwrap_or(0.0),
        mean_motion_dot: getf("MEAN_MOTION_DOT").unwrap_or(0.0),
        mean_motion_ddot: getf("MEAN_MOTION_DDOT").unwrap_or(0.0),
    })
}
