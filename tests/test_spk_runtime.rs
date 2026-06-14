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
