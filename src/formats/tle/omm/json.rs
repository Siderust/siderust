// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! JSON encoding for OMM, matching Celestrak's flat schema.

use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;
use serde::{Deserialize, Serialize};

use super::{format_epoch, parse_epoch, Omm};
use crate::formats::tle::TleError;
use crate::formats::tle::{Classification, SatelliteNumber};

/// Render an [`Omm`] (or a slice of them) to Celestrak-style JSON.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{omm::{Omm, json}, parse_3le};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let omm = Omm::from_tle(&parse_3le("ISS (ZARYA)", l1, l2).unwrap());
/// let s = json::write(&omm).unwrap();
/// assert!(s.contains("\"NORAD_CAT_ID\":25544"));
/// ```
pub fn write(omm: &Omm) -> Result<String, TleError> {
    let dto = OmmDto::from_omm(omm)?;
    serde_json::to_string(&dto).map_err(|e| TleError::OmmJson(e.to_string()))
}

/// Render a slice of [`Omm`] as a JSON array.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{omm::{Omm, json}, parse_3le};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let omm = Omm::from_tle(&parse_3le("ISS (ZARYA)", l1, l2).unwrap());
/// let arr = json::write_many(&[omm]).unwrap();
/// assert!(arr.starts_with('['));
/// ```
pub fn write_many(items: &[Omm]) -> Result<String, TleError> {
    let dtos: Vec<OmmDto> = items
        .iter()
        .map(OmmDto::from_omm)
        .collect::<Result<_, _>>()?;
    serde_json::to_string(&dtos).map_err(|e| TleError::OmmJson(e.to_string()))
}

/// Parse a single [`Omm`] from a JSON object.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::omm::json;
/// let s = r#"{"OBJECT_NAME":"X","OBJECT_ID":"1998-067A",
///   "EPOCH":"2008-09-20T12:25:40.104192","MEAN_MOTION":15.72125391,
///   "ECCENTRICITY":0.0006703,"INCLINATION":51.6416,
///   "RA_OF_ASC_NODE":247.4627,"ARG_OF_PERICENTER":130.5360,
///   "MEAN_ANOMALY":325.0288,"EPHEMERIS_TYPE":0,
///   "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":25544,
///   "ELEMENT_SET_NO":292,"REV_AT_EPOCH":56353,
///   "BSTAR":-1.1606e-5,"MEAN_MOTION_DOT":-2.182e-5,"MEAN_MOTION_DDOT":0.0}"#;
/// let omm = json::read(s).unwrap();
/// assert_eq!(omm.norad_id.0, 25544);
/// ```
pub fn read(input: &str) -> Result<Omm, TleError> {
    let dto: OmmDto = serde_json::from_str(input).map_err(|e| TleError::OmmJson(e.to_string()))?;
    dto.into_omm()
}

/// Parse a JSON array of OMMs.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::omm::json;
/// let arr = "[]";
/// let omms = json::read_many(arr).unwrap();
/// assert!(omms.is_empty());
/// ```
pub fn read_many(input: &str) -> Result<Vec<Omm>, TleError> {
    let dtos: Vec<OmmDto> =
        serde_json::from_str(input).map_err(|e| TleError::OmmJson(e.to_string()))?;
    dtos.into_iter().map(OmmDto::into_omm).collect()
}

/// Wire-shape DTO for the OMM JSON document. Field names match the
/// Celestrak GP catalog payload (uppercase, no transformation).
#[derive(Serialize, Deserialize)]
#[allow(non_snake_case)]
struct OmmDto {
    OBJECT_NAME: String,
    OBJECT_ID: String,
    EPOCH: String,
    MEAN_MOTION: f64,
    ECCENTRICITY: f64,
    INCLINATION: f64,
    RA_OF_ASC_NODE: f64,
    ARG_OF_PERICENTER: f64,
    MEAN_ANOMALY: f64,
    #[serde(default)]
    EPHEMERIS_TYPE: u8,
    #[serde(default = "default_classification")]
    CLASSIFICATION_TYPE: String,
    NORAD_CAT_ID: u32,
    #[serde(default)]
    ELEMENT_SET_NO: u16,
    #[serde(default)]
    REV_AT_EPOCH: u32,
    #[serde(default)]
    BSTAR: f64,
    #[serde(default)]
    MEAN_MOTION_DOT: f64,
    #[serde(default)]
    MEAN_MOTION_DDOT: f64,
    #[serde(default)]
    MEAN_ELEMENT_THEORY: Option<String>,
}

fn default_classification() -> String {
    "U".into()
}

impl OmmDto {
    fn from_omm(omm: &Omm) -> Result<Self, TleError> {
        Ok(Self {
            OBJECT_NAME: omm.object_name.clone(),
            OBJECT_ID: omm.object_id.clone(),
            EPOCH: format_epoch(omm.epoch)?,
            MEAN_MOTION: omm.mean_motion.value(),
            ECCENTRICITY: omm.eccentricity,
            INCLINATION: omm.inclination.value(),
            RA_OF_ASC_NODE: omm.ra_of_asc_node.value(),
            ARG_OF_PERICENTER: omm.arg_of_pericenter.value(),
            MEAN_ANOMALY: omm.mean_anomaly.value(),
            EPHEMERIS_TYPE: omm.ephemeris_type,
            CLASSIFICATION_TYPE: omm.classification.as_char().to_string(),
            NORAD_CAT_ID: omm.norad_id.0,
            ELEMENT_SET_NO: omm.element_set_no,
            REV_AT_EPOCH: omm.rev_at_epoch,
            BSTAR: omm.bstar,
            MEAN_MOTION_DOT: omm.mean_motion_dot,
            MEAN_MOTION_DDOT: omm.mean_motion_ddot,
            MEAN_ELEMENT_THEORY: Some("SGP4".to_string()),
        })
    }

    fn into_omm(self) -> Result<Omm, TleError> {
        if let Some(t) = &self.MEAN_ELEMENT_THEORY {
            if !t.eq_ignore_ascii_case("SGP4") {
                return Err(TleError::OmmUnsupportedTheory(t.clone()));
            }
        }
        let cls_char = self.CLASSIFICATION_TYPE.chars().next().unwrap_or('U');
        Ok(Omm {
            object_name: self.OBJECT_NAME,
            object_id: self.OBJECT_ID,
            epoch: parse_epoch(&self.EPOCH)?,
            mean_motion: AngularRate::<Turn, Day>::new(self.MEAN_MOTION),
            eccentricity: self.ECCENTRICITY,
            inclination: Degrees::new(self.INCLINATION),
            ra_of_asc_node: Degrees::new(self.RA_OF_ASC_NODE),
            arg_of_pericenter: Degrees::new(self.ARG_OF_PERICENTER),
            mean_anomaly: Degrees::new(self.MEAN_ANOMALY),
            ephemeris_type: self.EPHEMERIS_TYPE,
            classification: Classification::from_char(cls_char)?,
            norad_id: SatelliteNumber(self.NORAD_CAT_ID),
            element_set_no: self.ELEMENT_SET_NO,
            rev_at_epoch: self.REV_AT_EPOCH,
            bstar: self.BSTAR,
            mean_motion_dot: self.MEAN_MOTION_DOT,
            mean_motion_ddot: self.MEAN_MOTION_DDOT,
        })
    }
}
