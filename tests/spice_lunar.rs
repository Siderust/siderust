// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon
//! Independent CSPICE DE440 lunar geometry regressions.
#![cfg(feature = "spice")]
use siderust::qtty::{DEG, M};
use siderust::spice::{LunarObserver, LunarTarget, SpiceContext};
use siderust::tempoch::{Time, TDB};

fn excerpt() -> SpiceContext {
    let mut ctx = SpiceContext::new();
    for (alias, bytes) in [
        (
            "de440s excerpt",
            include_bytes!("fixtures/lunar_de440/excerpt.bsp").as_slice(),
        ),
        (
            "DE440 BPC excerpt",
            include_bytes!("fixtures/lunar_de440/excerpt.bpc").as_slice(),
        ),
        (
            "moon_de440_220930.tf",
            include_bytes!("fixtures/lunar_de440/moon_de440_220930.tf").as_slice(),
        ),
        (
            "naif0012.tls",
            include_bytes!("fixtures/lunar_de440/naif0012.tls").as_slice(),
        ),
    ] {
        ctx.load_bytes(alias, bytes.to_vec()).unwrap();
    }
    ctx.load_bytes(
        "synthetic large rotation",
        include_bytes!("fixtures/lunar_de440/test_rotation.tf").to_vec(),
    )
    .unwrap();
    ctx
}

fn compare(ctx: &SpiceContext, reference: &str) {
    let mut maxima = [0.0_f64; 5];
    for (line, text) in reference.lines().enumerate().skip(1) {
        let v: Vec<f64> = text.split(',').map(|x| x.parse().unwrap()).collect();
        let frame = ctx
            .frame_registry()
            .frame_by_id(v[5] as i32)
            .unwrap()
            .name
            .clone();
        let observer =
            LunarObserver::radial(ctx, &frame, v[1] * DEG, v[2] * DEG, v[3] * 1000.0 * M).unwrap();
        let target = if v[4] == 10.0 {
            LunarTarget::Sun
        } else {
            LunarTarget::Earth
        };
        let obs = observer.observe(target, Time::<TDB>::new(v[0])).unwrap();
        let m = ctx.rotation_naif(1, v[5] as i32, v[0]).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let err = (m[i][j] - v[6 + 3 * i + j]).abs();
                maxima[0] = maxima[0].max(err);
                assert!(err < 5e-14, "line {line} matrix ({i},{j}): {err}");
                let dot: f64 = (0..3).map(|k| m[i][k] * m[j][k]).sum();
                assert!((dot - if i == j { 1.0 } else { 0.0 }).abs() < 2e-14);
            }
        }
        for i in 0..3 {
            let err = (obs.target_fixed[i].value() - v[15 + i]).abs();
            maxima[1] = maxima[1].max(err);
            assert!(err < 2e-7, "line {line} target: {err}");
            let err = (obs.relative_fixed[i].value() - v[18 + i]).abs();
            maxima[2] = maxima[2].max(err);
            assert!(err < 2e-7, "line {line} relative: {err}");
        }
        let alt = (obs.horizontal.alt().value() - v[21]).abs();
        let az = ((obs.horizontal.az().value() - v[22] + 180.0).rem_euclid(360.0) - 180.0).abs();
        maxima[3] = maxima[3].max(alt.max(az));
        assert!(
            alt < 1e-10 && az < 1e-10,
            "line {line} angles alt={alt} az={az}"
        );
        let range = (obs.horizontal.range().value() - v[23]).abs();
        maxima[4] = maxima[4].max(range);
        assert!(range < 2e-7, "line {line} range {range}");
        assert_eq!(obs.above_horizon(), v[24] == 1.0, "line {line}");
        assert_eq!(obs.frame_id, v[5] as i32);
        assert!(!obs.kernels.is_empty());
    }
    eprintln!(
        "CSPICE maxima: matrix={}, target km={}, relative km={}, angle deg={}, range km={}",
        maxima[0], maxima[1], maxima[2], maxima[3], maxima[4]
    );
}

#[test]
fn cspice_de440_excerpt_parity() {
    compare(
        &excerpt(),
        include_str!("fixtures/lunar_de440/reference.csv"),
    );
}

#[test]
#[ignore = "Set SIDERUST_LUNAR_KERNEL_DIR to the five original NAIF kernels"]
fn cspice_external_kernel_parity() {
    let directory = std::env::var("SIDERUST_LUNAR_KERNEL_DIR").unwrap();
    let mut ctx = SpiceContext::new();
    for f in [
        "de440s.bsp",
        "moon_pa_de440_200625.bpc",
        "moon_de440_220930.tf",
        "pck00011.tpc",
        "naif0012.tls",
    ] {
        ctx.load(std::path::Path::new(&directory).join(f)).unwrap();
    }
    compare(
        &ctx,
        include_str!("fixtures/lunar_de440/external_reference.csv"),
    );
    for text in include_str!("fixtures/lunar_de440/coverage_rotations.csv").lines() {
        let v: Vec<f64> = text.split(',').map(|s| s.parse().unwrap()).collect();
        let m = ctx.rotation_naif(1, 31008, v[0]).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((m[i][j] - v[1 + i * 3 + j]).abs() < 5e-12);
            }
        }
    }
}

#[test]
fn frame_inverses_continuity_and_explicit_errors() {
    let ctx = excerpt();
    let start = 819720000.0;
    for et in [start, start + 691200.0, start + 2.0 * 691200.0, 822484800.0] {
        let forward = ctx.rotation_naif(1, 31011, et).unwrap();
        let back = ctx.rotation_naif(31011, 1, et).unwrap();
        let det = forward[0][0] * (forward[1][1] * forward[2][2] - forward[1][2] * forward[2][1])
            - forward[0][1] * (forward[1][0] * forward[2][2] - forward[1][2] * forward[2][0])
            + forward[0][2] * (forward[1][0] * forward[2][1] - forward[1][1] * forward[2][0]);
        assert!((det - 1.0).abs() < 2e-14);
        for i in 0..3 {
            for j in 0..3 {
                assert!((forward[i][j] - back[j][i]).abs() < 1e-14);
            }
        }
        if et > start && et < 822484800.0 {
            let a = ctx.rotation_naif(1, 31011, et - 0.001).unwrap();
            let b = ctx.rotation_naif(1, 31011, et + 0.001).unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    assert!((a[i][j] - b[i][j]).abs() < 1e-8);
                }
            }
        }
    }
    assert!(ctx.rotation_naif(1, 31011, start - 1.0).is_err());
    assert!(ctx.rotation_naif(1, 31011, 822484800.0 + 1.0).is_err());
    assert!(ctx.rotation_naif(1, 31011, f64::INFINITY).is_err());
    assert!(LunarObserver::radial(&ctx, "J2000", 0.0 * DEG, 0.0 * DEG, 1737400.0 * M).is_err());
    assert!(LunarObserver::radial(&ctx, "MOON_ME", 0.0 * DEG, 0.0 * DEG, -1.0 * M).is_err());
    assert!(LunarObserver::radial(&ctx, "IAU_MOON", 0.0 * DEG, 0.0 * DEG, 1737400.0 * M).is_err());
    assert!(SpiceContext::new().rotation_naif(1, 31011, start).is_err());
    let mut cycle = SpiceContext::new();
    cycle.load_bytes("cycle",b"\\begindata\nFRAME_77_NAME = 'CYCLE'\nFRAME_77_CLASS = 4\nFRAME_77_CLASS_ID = 77\nFRAME_77_CENTER = 301\nTKFRAME_77_SPEC = 'MATRIX'\nTKFRAME_77_RELATIVE = 'CYCLE'\nTKFRAME_77_MATRIX = (1 0 0 0 1 0 0 0 1)\n".to_vec()).unwrap();
    assert!(cycle.rotation_naif(1, 77, start).is_err());
    assert!(LunarObserver::radial(&cycle, "CYCLE", 0.0 * DEG, 0.0 * DEG, 1737400.0 * M).is_err());
}

#[test]
fn spk_stack_respects_per_target_priority_across_files() {
    let mut ctx = excerpt();
    let et = 819720000.0;
    let baseline = ctx.state_naif(10, 301, et).unwrap();
    let mut b = include_bytes!("fixtures/lunar_de440/excerpt.bsp").to_vec();
    let daf = siderust::formats::spice::daf::Daf::parse(&b).unwrap();
    let sun = &daf.summaries[0];
    assert_eq!(sun.target_id, 10);
    // Leave only Sun in this file, requiring older kernel edges for the Moon.
    let record = i32::from_le_bytes(b[76..80].try_into().unwrap()) as usize;
    let off = (record - 1) * 1024;
    b[off + 16..off + 24].copy_from_slice(&1.0_f64.to_le_bytes());
    let word = sun.start_word + 2;
    let off = (word - 1) * 8;
    let value = f64::from_le_bytes(b[off..off + 8].try_into().unwrap());
    b[off..off + 8].copy_from_slice(&(value + 1.0).to_le_bytes());
    ctx.load_bytes("Sun override", b).unwrap();
    let changed = ctx.state_naif(10, 301, et).unwrap();
    assert!((changed[0] - baseline[0] - 1.0).abs() < 1e-7);
    for i in 1..6 {
        assert_eq!(changed[i], baseline[i]);
    }
}

#[test]
fn bpc_priority_never_falls_back_on_unsupported_or_uncovered_class() {
    let mut ctx = excerpt();
    let mut b = include_bytes!("fixtures/lunar_de440/excerpt.bpc").to_vec();
    let rec = i32::from_le_bytes(b[76..80].try_into().unwrap()) as usize;
    let off = (rec - 1) * 1024 + 24;
    // High-priority file has a narrower interval and unsupported Type 20.
    b[off..off + 8].copy_from_slice(&819721000.0_f64.to_le_bytes());
    b[off + 24..off + 28].copy_from_slice(&20i32.to_le_bytes());
    ctx.load_bytes("unsupported overriding BPC", b.clone())
        .unwrap();
    assert!(ctx.rotation_naif(1, 31011, 819721001.0).is_err());
    assert!(ctx.rotation_naif(1, 31011, 819720000.0).is_ok());
    // Also verify later segment wins within one file.
    let mut b = include_bytes!("fixtures/lunar_de440/excerpt.bpc").to_vec();
    let desc = b[off..off + 40].to_vec();
    b[off + 40..off + 80].copy_from_slice(&desc);
    b[off - 8..off].copy_from_slice(&2.0_f64.to_le_bytes());
    b[off + 40 + 24..off + 40 + 28].copy_from_slice(&20i32.to_le_bytes());
    let mut ctx = excerpt();
    ctx.load_bytes("two segments", b).unwrap();
    assert!(ctx.rotation_naif(1, 31011, 819720000.0).is_err());
}

#[test]
fn typed_utc_tt_tdb_with_historical_leaps() {
    use siderust::tempoch::{TimeContext, TT, UTC};
    let context = TimeContext::new();
    for line in include_str!("fixtures/lunar_de440/time_reference.csv").lines() {
        let cols: Vec<_> = line.split(',').collect();
        let utc = Time::<UTC>::try_from_chrono_with(
            cols[0].parse::<chrono::DateTime<chrono::Utc>>().unwrap(),
            &context,
        )
        .unwrap();
        let tt = utc.to::<TT>().to_j2000s().raw().value();
        let tdb = utc.to::<TDB>().to_j2000s().raw().value();
        let reference_tt = cols[1].parse::<f64>().unwrap();
        let reference_et = cols[2].parse::<f64>().unwrap();
        assert!(
            (tt - reference_tt).abs() < 1e-6,
            "{} TT residual {}",
            cols[0],
            tt - reference_tt
        );
        // tempoch and LSK DELTET have distinct TDB approximations; quantify it
        // separately from the exact-ET geometry comparisons above.
        assert!(
            (tdb - reference_et).abs() < 5e-5,
            "{} TDB residual {}",
            cols[0],
            tdb - reference_et
        );
    }
}

#[test]
fn malformed_spk_descriptors_and_directories_are_rejected() {
    let good = include_bytes!("fixtures/lunar_de440/excerpt.bsp").to_vec();
    let daf = siderust::formats::spice::daf::Daf::parse(&good).unwrap();
    let rec = i32::from_le_bytes(good[76..80].try_into().unwrap()) as usize;
    let desc = (rec - 1) * 1024 + 24;
    for offset in [desc + 32, desc + 36] {
        for value in [-1, i32::MAX] {
            let mut b = good.clone();
            b[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
            assert!(SpiceContext::new().load_bytes("bad SPK", b).is_err());
        }
    }
    let end = daf.summaries[0].end_word;
    for (word, value) in [
        (end, f64::NAN),
        (end, 1.5),
        (end - 1, 8.5),
        (end - 1, 1e300),
        (daf.summaries[0].start_word + 2, f64::NAN),
    ] {
        let mut b = good.clone();
        let off = (word - 1) * 8;
        b[off..off + 8].copy_from_slice(&value.to_le_bytes());
        assert!(SpiceContext::new().load_bytes("bad SPK", b).is_err());
    }
}
