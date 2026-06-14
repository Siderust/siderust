// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Comprehensive unit tests for `siderust-tle`.
//!
//! Drives line-coverage of every parser branch (TLE column-extractor,
//! Alpha-5, checksum, OMM KVN/XML/JSON) and round-trip parity of the
//! three OMM encodings.

#[cfg(feature = "serde")]
use super::omm::json;
use super::omm::{kvn, xml, Omm};
use super::parse::{compute_checksum, parse_assumed_decimal_exponent};
use super::*;

const ISS_NAME: &str = "ISS (ZARYA)";
const ISS_L1: &str = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
const ISS_L2: &str = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

// A second canonical TLE (Vanguard 1, lowest NORAD id), for variety.
const VAN_L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const VAN_L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

#[test]
fn checksum_of_canonical_lines_validates() {
    validate_tle_checksum(ISS_L1).unwrap();
    validate_tle_checksum(ISS_L2).unwrap();
    validate_tle_checksum(VAN_L1).unwrap();
    validate_tle_checksum(VAN_L2).unwrap();
}

#[test]
fn parse_iss_returns_expected_fields() {
    let tle = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
    assert_eq!(tle.norad_id, SatelliteNumber(25544));
    assert_eq!(tle.classification, Classification::Unclassified);
    assert_eq!(tle.international_designator.0, "98067A");
    assert!((tle.inclination.value() - 51.6416).abs() < 1e-6);
    assert!((tle.raan.value() - 247.4627).abs() < 1e-6);
    assert!((tle.eccentricity - 0.0006703).abs() < 1e-9);
    assert!((tle.argument_of_perigee.value() - 130.5360).abs() < 1e-6);
    assert!((tle.mean_anomaly.value() - 325.0288).abs() < 1e-6);
    assert!((tle.mean_motion.value() - 15.72125391).abs() < 1e-8);
    assert!((tle.bstar - -0.11606e-4).abs() < 1e-12);
    assert_eq!(tle.element_set_number, 292);
    assert_eq!(tle.revolution_number_at_epoch, 56353);
}

#[test]
fn parse_3le_strips_zero_prefix() {
    let with_prefix = parse_3le("0 ISS (ZARYA)", ISS_L1, ISS_L2).unwrap();
    assert_eq!(with_prefix.name.as_deref(), Some("ISS (ZARYA)"));
    let empty = parse_3le("   ", ISS_L1, ISS_L2).unwrap();
    assert!(empty.name.is_none());
}

#[test]
fn invalid_checksum_is_rejected() {
    let mut bad = ISS_L1.to_string();
    let n = bad.len();
    bad.replace_range(n - 1..n, "0");
    match validate_tle_checksum(&bad) {
        Err(TleError::BadChecksum { line: 1, .. }) => {}
        other => panic!("expected BadChecksum on line 1, got {other:?}"),
    }
}

#[test]
fn non_digit_checksum_is_rejected() {
    let mut bad = ISS_L1.to_string();
    let n = bad.len();
    bad.replace_range(n - 1..n, "X");
    assert!(matches!(
        validate_tle_checksum(&bad),
        Err(TleError::BadChecksum { .. })
    ));
}

#[test]
fn bad_length_is_rejected() {
    let short = &ISS_L1[..68];
    assert!(matches!(
        validate_tle_checksum(short),
        Err(TleError::BadLength { .. })
    ));
    let short = &ISS_L1[..50];
    assert!(matches!(
        parse_tle(short, ISS_L2),
        Err(TleError::BadLength { line: 1, .. })
    ));
    assert!(matches!(
        parse_tle(ISS_L1, &ISS_L2[..50]),
        Err(TleError::BadLength { line: 2, .. })
    ));
}

#[test]
fn bad_leading_char_is_rejected() {
    // line1 with leading '9' must produce a checksum length error first because we
    // reject on length parity above; craft a 69-char line that just changes col 0.
    let mut bad1 = ISS_L1.to_string();
    bad1.replace_range(0..1, "9");
    // Also fix checksum so we hit BadLeadingChar specifically.
    let cs = compute_tle_checksum(&bad1[..68]);
    bad1.replace_range(68..69, &cs.to_string());
    assert!(matches!(
        parse_tle(&bad1, ISS_L2),
        Err(TleError::BadLeadingChar { line: 1, .. })
    ));

    let mut bad2 = ISS_L2.to_string();
    bad2.replace_range(0..1, "9");
    let cs = compute_tle_checksum(&bad2[..68]);
    bad2.replace_range(68..69, &cs.to_string());
    assert!(matches!(
        parse_tle(ISS_L1, &bad2),
        Err(TleError::BadLeadingChar { line: 2, .. })
    ));
}

#[test]
fn alpha5_decoding_full_table() {
    // Regular numeric.
    assert_eq!(
        SatelliteNumber::parse("25544").unwrap(),
        SatelliteNumber(25_544)
    );
    // Each Alpha-5 segment boundary.
    assert_eq!(
        SatelliteNumber::parse("A0000").unwrap(),
        SatelliteNumber(100_000)
    );
    assert_eq!(
        SatelliteNumber::parse("H9999").unwrap(),
        SatelliteNumber(179_999)
    );
    assert_eq!(
        SatelliteNumber::parse("J0000").unwrap(),
        SatelliteNumber(180_000)
    );
    assert_eq!(
        SatelliteNumber::parse("N9999").unwrap(),
        SatelliteNumber(229_999)
    );
    assert_eq!(
        SatelliteNumber::parse("P0000").unwrap(),
        SatelliteNumber(230_000)
    );
    assert_eq!(
        SatelliteNumber::parse("T0001").unwrap(),
        SatelliteNumber(270_001)
    );
    assert_eq!(
        SatelliteNumber::parse("Z9999").unwrap(),
        SatelliteNumber(339_999)
    );

    // Forbidden letters and shape errors.
    assert!(matches!(
        SatelliteNumber::parse("I0001"),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse("O0001"),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse(""),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse("A001"),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse("A00X1"),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse("aaaaa"),
        Err(TleError::InvalidAlpha5 { .. })
    ));
    assert!(matches!(
        SatelliteNumber::parse("9".repeat(6).as_str()),
        Err(TleError::InvalidAlpha5 { .. })
    ));
}

#[test]
fn alpha5_format_round_trip() {
    for n in [
        0u32, 5, 25_544, 99_999, 100_000, 179_999, 180_000, 270_001, 339_999,
    ] {
        let s = SatelliteNumber(n).format_alpha5().unwrap();
        assert_eq!(
            SatelliteNumber::parse(&s).unwrap(),
            SatelliteNumber(n),
            "round-trip {n}"
        );
    }
    assert!(SatelliteNumber(340_000).format_alpha5().is_err());
}

#[test]
fn classification_round_trip() {
    for &c in &['U', 'C', 'S'] {
        let parsed = Classification::from_char(c).unwrap();
        assert_eq!(parsed.as_char(), c);
    }
    assert!(matches!(
        Classification::from_char('Z'),
        Err(TleError::InvalidClassification { .. })
    ));
}

#[test]
fn assumed_decimal_with_exponent_branches() {
    assert!((parse_assumed_decimal_exponent(" 12345-3", "x").unwrap() - 0.12345e-3).abs() < 1e-12);
    assert!((parse_assumed_decimal_exponent("-11606-4", "x").unwrap() - -0.11606e-4).abs() < 1e-15);
    assert!((parse_assumed_decimal_exponent("+12345+2", "x").unwrap() - 12.345).abs() < 1e-12);
    assert_eq!(parse_assumed_decimal_exponent("00000-0", "x").unwrap(), 0.0);
    assert_eq!(
        parse_assumed_decimal_exponent("        ", "x").unwrap(),
        0.0
    );

    // Failure modes.
    assert!(parse_assumed_decimal_exponent("-", "x").is_err());
    assert!(parse_assumed_decimal_exponent("ab", "x").is_err());
    assert!(parse_assumed_decimal_exponent("abcde-3", "x").is_err());
    assert!(parse_assumed_decimal_exponent("1234567", "x").is_err());
}

#[test]
fn checksum_helper_matches_bytes() {
    // Computed checksum equals the digit baked into the canonical fixture.
    assert_eq!(compute_checksum(&ISS_L1[..68]), b'7' - b'0');
    assert_eq!(compute_checksum(&ISS_L2[..68]), b'7' - b'0');
    assert_eq!(compute_tle_checksum(&VAN_L1[..68]), b'3' - b'0');
}

#[test]
fn epoch_year_expansion_and_invalid_epoch() {
    use super::parse::expand_two_digit_year;
    assert_eq!(expand_two_digit_year(57), 1957);
    assert_eq!(expand_two_digit_year(99), 1999);
    assert_eq!(expand_two_digit_year(00), 2000);
    assert_eq!(expand_two_digit_year(56), 2056);

    use super::parse::epoch_from_year_doy;
    assert!(matches!(
        epoch_from_year_doy(2024, 0.0),
        Err(TleError::InvalidEpoch { .. })
    ));
    assert!(matches!(
        epoch_from_year_doy(2024, 367.5),
        Err(TleError::InvalidEpoch { .. })
    ));
    assert!(matches!(
        epoch_from_year_doy(2023, 366.0),
        Err(TleError::InvalidEpoch { .. })
    ));
    assert!(matches!(
        epoch_from_year_doy(2024, f64::NAN),
        Err(TleError::InvalidEpoch { .. })
    ));
    // Valid: 2024 was a leap year.
    epoch_from_year_doy(2024, 366.5).unwrap();
}

#[test]
fn mismatched_norad_id_is_rejected() {
    let l2_other = "2 25545  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563530";
    let computed = compute_tle_checksum(&l2_other[..68]);
    let mut forged = l2_other.to_string();
    forged.replace_range(68..69, &computed.to_string());
    assert!(matches!(
        parse_tle(ISS_L1, &forged),
        Err(TleError::MismatchedSatelliteNumber { .. })
    ));
}

#[test]
fn invalid_classification_field_is_rejected() {
    let mut bad = ISS_L1.to_string();
    bad.replace_range(7..8, "Z");
    let cs = compute_tle_checksum(&bad[..68]);
    bad.replace_range(68..69, &cs.to_string());
    assert!(matches!(
        parse_tle(&bad, ISS_L2),
        Err(TleError::InvalidClassification { .. })
    ));
}

#[test]
fn invalid_eccentricity_field_is_rejected() {
    // Replace the eccentricity columns 26..33 with non-digits.
    let mut bad = ISS_L2.to_string();
    bad.replace_range(26..33, "ABCDEFG");
    let cs = compute_tle_checksum(&bad[..68]);
    bad.replace_range(68..69, &cs.to_string());
    assert!(matches!(
        parse_tle(ISS_L1, &bad),
        Err(TleError::InvalidNumber {
            field: "eccentricity",
            ..
        })
    ));
}

#[test]
fn builder_round_trip_to_format_tle() {
    use super::builder::format_tle;
    let tle = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    let (l1, l2) = format_tle(&tle).unwrap();
    assert_eq!(l1.len(), 69);
    assert_eq!(l2.len(), 69);
    // Re-parse and compare deeply (the byte-for-byte encoding can vary
    // in zero-padding, but the parsed values must be identical).
    let tle2 = parse_tle(&l1, &l2).unwrap();
    assert_eq!(tle.norad_id, tle2.norad_id);
    assert_eq!(tle.classification, tle2.classification);
    assert_eq!(tle.international_designator, tle2.international_designator);
    assert!((tle.inclination.value() - tle2.inclination.value()).abs() < 1e-4);
    assert!((tle.mean_motion.value() - tle2.mean_motion.value()).abs() < 1e-8);
    assert!((tle.bstar - tle2.bstar).abs() < 1e-9);
    assert_eq!(tle.element_set_number, tle2.element_set_number);
    assert_eq!(
        tle.revolution_number_at_epoch,
        tle2.revolution_number_at_epoch
    );
}

#[test]
fn builder_missing_field_errors() {
    let err = TleBuilder::new().build().unwrap_err();
    assert!(matches!(err, TleError::BuilderMissingField(_)));
}

#[test]
fn omm_kvn_round_trip() {
    let tle = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    let omm = Omm::from_tle(&tle);
    let text = kvn::write(&omm).unwrap();
    let omm2 = kvn::read(&text).unwrap();
    assert_eq!(omm2.norad_id, omm.norad_id);
    assert_eq!(omm2.object_id, omm.object_id);
    assert!((omm2.inclination.value() - omm.inclination.value()).abs() < 1e-6);
    assert!((omm2.bstar - omm.bstar).abs() < 1e-12);
    assert_eq!(omm2.classification, omm.classification);
}

#[test]
fn omm_xml_round_trip() {
    let tle = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    let omm = Omm::from_tle(&tle);
    let s = xml::write(&omm).unwrap();
    let omm2 = xml::read(&s).unwrap();
    assert_eq!(omm2.norad_id, omm.norad_id);
    assert!((omm2.eccentricity - omm.eccentricity).abs() < 1e-9);
}

#[cfg(feature = "serde")]
#[test]
fn omm_json_round_trip_single_and_array() {
    let tle = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    let omm = Omm::from_tle(&tle);
    let s = json::write(&omm).unwrap();
    let omm2 = json::read(&s).unwrap();
    assert_eq!(omm2.norad_id, omm.norad_id);

    let arr = json::write_many(&[omm.clone(), omm.clone()]).unwrap();
    let parsed = json::read_many(&arr).unwrap();
    assert_eq!(parsed.len(), 2);
}

#[test]
fn omm_kvn_rejects_unsupported_theory() {
    let bad = "MEAN_ELEMENT_THEORY = SDP4\nOBJECT_NAME = X\n";
    assert!(matches!(
        kvn::read(bad),
        Err(TleError::OmmUnsupportedTheory(_))
    ));
}

#[test]
fn omm_kvn_missing_field() {
    let bad = "MEAN_ELEMENT_THEORY = SGP4\n";
    assert!(matches!(kvn::read(bad), Err(TleError::OmmMissingField(_))));
}

#[test]
fn omm_kvn_invalid_number() {
    // All required fields present; INCLINATION malformed.
    let bad = "OBJECT_NAME = X\nOBJECT_ID = 1998-067A\n\
               EPOCH = 2008-09-20T12:25:40.000\nMEAN_MOTION = 15.7\n\
               ECCENTRICITY = 0.001\nINCLINATION = NOT_A_NUMBER\n\
               RA_OF_ASC_NODE = 0\nARG_OF_PERICENTER = 0\nMEAN_ANOMALY = 0\n\
               NORAD_CAT_ID = 25544\n";
    assert!(matches!(kvn::read(bad), Err(TleError::OmmMalformed { .. })));
}

#[test]
fn omm_xml_missing_field() {
    let bad = "<omm><body><segment><data><meanElements></meanElements>\
               </data></segment></body></omm>";
    assert!(matches!(xml::read(bad), Err(TleError::OmmMissingField(_))));
}

#[cfg(feature = "serde")]
#[test]
fn omm_json_invalid_returns_json_error() {
    assert!(matches!(
        json::read("{ malformed"),
        Err(TleError::OmmJson(_))
    ));
}

#[test]
fn omm_round_trip_to_tle() {
    let tle1 = parse_3le(ISS_NAME, ISS_L1, ISS_L2).unwrap();
    let omm = Omm::from_tle(&tle1);
    let tle2 = omm.to_tle();
    assert_eq!(tle1.norad_id, tle2.norad_id);
    assert_eq!(tle1.international_designator, tle2.international_designator);
    assert!((tle1.inclination.value() - tle2.inclination.value()).abs() < 1e-9);
    assert!((tle1.mean_motion.value() - tle2.mean_motion.value()).abs() < 1e-12);
}

#[test]
fn omm_intl_designator_helpers() {
    use super::omm::{contract_intl_designator, expand_intl_designator};
    assert_eq!(expand_intl_designator("98067A"), "1998-067A");
    assert_eq!(expand_intl_designator("00179B"), "2000-179B");
    assert_eq!(contract_intl_designator("1998-067A"), "98067A");
    assert_eq!(contract_intl_designator("2000-179B"), "00179B");
    // Unparseable inputs round-trip unchanged.
    assert_eq!(expand_intl_designator("FOO"), "FOO");
    assert_eq!(contract_intl_designator("NOPE"), "NOPE");
}

#[test]
fn omm_xml_xss_escape_round_trip() {
    let tle = parse_3le("WEIRD<&\"NAME>", ISS_L1, ISS_L2).unwrap();
    let mut omm = Omm::from_tle(&tle);
    omm.object_name = "WEIRD<&\"NAME>".to_string();
    let s = xml::write(&omm).unwrap();
    let omm2 = xml::read(&s).unwrap();
    assert_eq!(omm2.object_name, omm.object_name);
}

#[test]
fn fixtures_round_trip() {
    let txt = include_str!("../test-data/tle/iss_zarya.3le");
    let mut iter = txt.lines();
    let name = iter.next().unwrap();
    let l1 = iter.next().unwrap();
    let l2 = iter.next().unwrap();
    let tle = parse_3le(name, l1, l2).unwrap();
    assert_eq!(tle.norad_id, SatelliteNumber(25544));

    let kvn_txt = include_str!("../test-data/omm/iss_zarya.kvn");
    let omm = kvn::read(kvn_txt).unwrap();
    assert_eq!(omm.norad_id, SatelliteNumber(25544));

    let xml_txt = include_str!("../test-data/omm/iss_zarya.xml");
    let omm_x = xml::read(xml_txt).unwrap();
    assert_eq!(omm_x.norad_id, SatelliteNumber(25544));
}

#[cfg(feature = "serde")]
#[test]
fn fixtures_round_trip_json() {
    let json_txt = include_str!("../test-data/omm/iss_zarya.json");
    let omm_j = json::read(json_txt).unwrap();
    assert_eq!(omm_j.norad_id, SatelliteNumber(25544));

    let arr_txt = include_str!("../test-data/omm/celestrak_sample.json");
    let many = json::read_many(arr_txt).unwrap();
    assert!(many.len() >= 2);
}
