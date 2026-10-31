// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime SPK / BSP integration tests (synthetic kernel + optional real BSP).

use siderust::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::formats::spice::spk::{self, IndexedSegmentData, SegmentData};
use siderust::formats::spice::SpkKernelSet;
use siderust::time::JulianDate;

const SECONDS_PER_DAY: f64 = siderust::qtty::time::SECONDS_PER_DAY;
const JD_J2000: f64 = tempoch::J2000_JD_TT_DAY.value();

fn segment_constant(x_km: f64, init: f64, intlen: f64) -> SegmentData {
    SegmentData {
        data_type: 2,
        init,
        intlen,
        rsize: 8,
        ncoeff: 2,
        n_records: 1,
        records: vec![
            init + intlen / 2.0,
            intlen / 2.0,
            x_km,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    }
}

fn indexed(target: i32, center: i32, x_km: f64, start_et: f64, end_et: f64) -> IndexedSegmentData {
    let intlen = end_et - start_et;
    IndexedSegmentData {
        target_id: target,
        center_id: center,
        frame_id: 1,
        start_et,
        end_et,
        data: segment_constant(x_km, start_et, intlen),
    }
}

fn de_chain_indexed() -> Vec<IndexedSegmentData> {
    let intlen = 1000.0 * SECONDS_PER_DAY;
    vec![
        indexed(spk::SUN_TARGET, spk::SUN_CENTER, 1.0e8, 0.0, intlen),
        indexed(spk::EMB_TARGET, spk::EMB_CENTER, 1.5e8, 0.0, intlen),
        indexed(spk::MOON_TARGET, spk::MOON_CENTER, 3.84e5, 0.0, intlen),
    ]
}

#[test]
fn runtime_ephemeris_multi_segment_boundary() {
    let intlen = 500.0 * SECONDS_PER_DAY;
    let indexed = vec![
        indexed(spk::SUN_TARGET, spk::SUN_CENTER, 100.0, 0.0, intlen),
        indexed(
            spk::SUN_TARGET,
            spk::SUN_CENTER,
            200.0,
            intlen,
            2.0 * intlen,
        ),
        indexed(spk::EMB_TARGET, spk::EMB_CENTER, 1.5e8, 0.0, 2.0 * intlen),
        indexed(
            spk::MOON_TARGET,
            spk::MOON_CENTER,
            3.84e5,
            0.0,
            2.0 * intlen,
        ),
    ];
    let eph = RuntimeEphemeris::from_indexed_segments(indexed);
    let jd_first = JulianDate::new(JD_J2000 + 250.0);
    let jd_second = JulianDate::new(JD_J2000 + 750.0);
    let p1 = eph.sun_barycentric(jd_first);
    let p2 = eph.sun_barycentric(jd_second);
    let au = 149_597_870.700;
    assert!((p1.x().value() * au - 100.0).abs() < 1.0);
    assert!((p2.x().value() * au - 200.0).abs() < 1.0);
}

#[test]
fn spk_kernel_set_matches_runtime_earth_moon_chain() {
    let indexed = de_chain_indexed();
    let eph = RuntimeEphemeris::from_indexed_segments(indexed);
    let set = SpkKernelSet::from_indexed_segments(de_chain_indexed());
    let jd = JulianDate::new(JD_J2000 + 500.0);
    let moon_rt = eph.moon_geocentric(jd);
    let moon_set = set
        .try_geometric_state(301, 399, jd)
        .expect("moon geocentric from kernel set");
    assert!(
        (moon_rt.x().value() - moon_set.x().value()).abs() < 1.0,
        "runtime vs SpkKernelSet Moon X mismatch"
    );
}

#[test]
fn optional_real_bsp_smoke() {
    let Ok(path) = std::env::var("SIDERUST_BSP_PATH") else {
        eprintln!("Skipping real BSP test: set SIDERUST_BSP_PATH");
        return;
    };
    let eph = RuntimeEphemeris::from_bsp(&path).expect("load BSP");
    let jd = siderust::J2000;
    let earth = eph.earth_barycentric(jd);
    assert!(earth.x().is_finite());
    assert!(eph.moon_geocentric(jd).x().is_finite());
}

#[test]
fn pre_parsed_non_j2000_frames_are_rejected() {
    for frame in [0, 2, 17, 10020, -1] {
        let mut segments = de_chain_indexed();
        segments[0].frame_id = frame;
        let error = RuntimeEphemeris::try_from_indexed_segments(segments).unwrap_err();
        assert!(error.to_string().contains(&format!("frame {frame}")));
        let mut segments = de_chain_indexed();
        segments[0].frame_id = frame;
        assert!(SpkKernelSet::try_from_indexed_segments(segments).is_err());
    }
}

#[test]
fn byte_loaders_reject_non_j2000_frame_before_reading_coefficients() {
    // A DAF Sun summary in ECLIPJ2000 (17). No coefficients are needed:
    // rejection must happen before any attempt to interpret them as ICRF.
    let mut bytes = vec![0u8; 2048];
    bytes[0..8].copy_from_slice(b"DAF/SPK ");
    bytes[8..12].copy_from_slice(&2i32.to_le_bytes());
    bytes[12..16].copy_from_slice(&6i32.to_le_bytes());
    bytes[76..80].copy_from_slice(&2i32.to_le_bytes());
    let summary = &mut bytes[1024..];
    summary[16..24].copy_from_slice(&1f64.to_le_bytes());
    summary[32..40].copy_from_slice(&86400f64.to_le_bytes());
    summary[40..44].copy_from_slice(&10i32.to_le_bytes());
    summary[48..52].copy_from_slice(&17i32.to_le_bytes());
    summary[52..56].copy_from_slice(&2i32.to_le_bytes());
    summary[56..60].copy_from_slice(&1i32.to_le_bytes());
    summary[60..64].copy_from_slice(&100i32.to_le_bytes());
    for error in [
        spk::parse_indexed_segments(&bytes).unwrap_err().to_string(),
        spk::parse_bsp(&bytes).unwrap_err().to_string(),
        RuntimeEphemeris::from_bytes(&bytes)
            .unwrap_err()
            .to_string(),
    ] {
        assert!(error.contains("frame 17"), "{error}");
    }
}

/// Complete synthetic BSP with constant-position Type 2 records. Descriptor
/// type/frame may deliberately be incompatible; loaders must reject them.
fn bsp_with_summaries(summaries: &[(i32, i32, i32)]) -> Vec<u8> {
    let mut bytes = vec![0u8; 4096];
    bytes[0..8].copy_from_slice(b"DAF/SPK ");
    bytes[8..12].copy_from_slice(&2i32.to_le_bytes());
    bytes[12..16].copy_from_slice(&6i32.to_le_bytes());
    bytes[76..80].copy_from_slice(&2i32.to_le_bytes());
    bytes[1040..1048].copy_from_slice(&(summaries.len() as f64).to_le_bytes());
    for (i, &(target, frame, kind)) in summaries.iter().enumerate() {
        let start_word = 385 + i * 12;
        let offset = 1048 + i * 40;
        bytes[offset + 8..offset + 16].copy_from_slice(&86400f64.to_le_bytes());
        for (j, value) in [
            target,
            0,
            frame,
            kind,
            start_word as i32,
            (start_word + 11) as i32,
        ]
        .iter()
        .enumerate()
        {
            let offset = offset + 16 + j * 4;
            bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
        }
        let data = (start_word - 1) * 8;
        for (j, value) in [
            43200.,
            43200.,
            100. + i as f64,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            86400.,
            8.,
            1.,
        ]
        .iter()
        .enumerate()
        {
            bytes[data + j * 8..data + (j + 1) * 8].copy_from_slice(&f64::to_le_bytes(*value));
        }
    }
    bytes
}

#[test]
fn strict_byte_policy_rejects_unrelated_or_higher_priority_incompatible_segments() {
    let valid = bsp_with_summaries(&[(10, 1, 2)]);
    assert!(RuntimeEphemeris::from_bytes(&valid).is_ok());
    for (target, frame, kind, diagnostic) in [
        (999, 17, 2, "frame 17"),
        (10, 17, 2, "frame 17"),
        (999, 1, 5, "Type 5"),
        (10, 1, 5, "Type 5"),
    ] {
        let bytes = bsp_with_summaries(&[(10, 1, 2), (target, frame, kind)]);
        assert!(RuntimeEphemeris::from_bytes(&bytes)
            .unwrap_err()
            .to_string()
            .contains(diagnostic));
        let mut set = SpkKernelSet::default();
        set.push_bytes(&valid).unwrap();
        assert!(set
            .push_bytes(&bytes)
            .unwrap_err()
            .to_string()
            .contains(diagnostic));
        assert_eq!(set.segment_count(), 1, "failed load must be atomic");
        assert!(spk::parse_bsp(&bytes)
            .unwrap_err()
            .to_string()
            .contains(diagnostic));
    }
}

#[test]
fn parsed_constructors_use_the_same_type_and_body_chain_policy() {
    let mut segments = de_chain_indexed();
    segments[0].data.data_type = 5;
    assert!(RuntimeEphemeris::try_from_indexed_segments(segments)
        .unwrap_err()
        .to_string()
        .contains("Type 5"));
    let mut segments = de_chain_indexed();
    segments[0].data.data_type = 5;
    assert!(SpkKernelSet::try_from_indexed_segments(segments).is_err());
    let mut segments = de_chain_indexed();
    segments.push(indexed(301, 399, 20., 0., 86400.));
    assert!(RuntimeEphemeris::try_from_indexed_segments(segments)
        .unwrap_err()
        .to_string()
        .contains("got center 399"));
}

#[test]
fn later_segments_win_only_within_their_declared_coverage() {
    let day = SECONDS_PER_DAY;
    let segments = vec![
        indexed(10, 0, 100., 0., 1000. * day),
        indexed(10, 0, 900., 400. * day, 600. * day),
    ];
    let eph = RuntimeEphemeris::try_from_indexed_segments(segments).unwrap();
    let set = SpkKernelSet::try_from_indexed_segments(vec![
        indexed(10, 0, 100., 0., 1000. * day),
        indexed(10, 0, 900., 400. * day, 600. * day),
    ])
    .unwrap();
    for (days, km) in [(250., 100.), (500., 900.), (750., 100.)] {
        let jd = JulianDate::new(JD_J2000 + days);
        let p = eph
            .sun_barycentric(jd)
            .to_unit::<siderust::qtty::Kilometer>();
        let norm = (p.x().value().powi(2) + p.y().value().powi(2) + p.z().value().powi(2)).sqrt();
        assert!((norm - km).abs() < 1e-9);
        let state = set.try_geometric_state(10, 0, jd).unwrap();
        assert!((state.x().value() - km).abs() < 1e-9);
    }
}

#[test]
fn broken_selected_coefficients_do_not_expose_older_segments() {
    let day = SECONDS_PER_DAY;
    let make = || {
        let mut newer = indexed(10, 0, 900., 0., 100. * day);
        newer.end_et = 1000. * day; // priority covers query, coefficients do not
        vec![indexed(10, 0, 100., 0., 1000. * day), newer]
    };
    let jd = JulianDate::new(JD_J2000 + 500.);
    let eph = RuntimeEphemeris::try_from_indexed_segments(make()).unwrap();
    assert!(eph.try_sun_barycentric(jd).is_err());
    let set = SpkKernelSet::try_from_indexed_segments(make()).unwrap();
    assert!(set.try_geometric_state(10, 0, jd).is_err());
}

#[test]
fn selected_target_chain_never_falls_back_to_an_older_center() {
    use siderust::formats::spice::SpkKernelError;
    let day = SECONDS_PER_DAY;
    let set = SpkKernelSet::try_from_indexed_segments(vec![
        indexed(10, 0, 100., 0., 1000. * day),
        indexed(10, 777, 900., 0., 1000. * day),
    ])
    .unwrap();
    assert!(matches!(
        set.try_geometric_state(10, 0, JulianDate::new(JD_J2000 + 500.)),
        Err(SpkKernelError::MissingSegment { target_id: 777 })
    ));
}
