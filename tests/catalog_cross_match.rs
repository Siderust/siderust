// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Public API tests for stellar catalog cross-matching.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use siderust::catalogs::{
    CatalogRecord, CrossMatchError, CrossMatchMode, CrossMatchOptions, CrossMatchSide,
    LargeStarCatalog, MissingProperMotion,
};
use siderust::qtty::{Arcseconds, JulianYears, Radian, Radians};

fn record(source_id: u64, ra: f64, dec: f64, epoch_jyr: f64) -> CatalogRecord {
    CatalogRecord {
        source_id,
        ra: Radians::new(ra),
        dec: Radians::new(dec),
        epoch_jyr,
        pm_ra_cosdec: None,
        pm_dec: None,
        parallax: None,
        radial_velocity: None,
        g_mag: None,
        bp_mag: None,
        rp_mag: None,
        quality_ok: true,
    }
}

fn catalog(rows: impl IntoIterator<Item = CatalogRecord>) -> LargeStarCatalog {
    let mut catalog = LargeStarCatalog::new();
    catalog.extend(rows);
    catalog
}

fn options(radius: Radians) -> CrossMatchOptions {
    CrossMatchOptions::new(JulianYears::new(2016.0), radius)
}

#[test]
fn zero_radius_matches_exact_directions_across_wrap_and_at_both_poles() {
    let left = catalog([
        record(10, TAU, 0.2, 2016.0),
        record(11, 1.0, FRAC_PI_2, 2016.0),
        record(12, 2.0, -FRAC_PI_2, 2016.0),
        record(13, 0.0, 0.0, 2016.0),
    ]);
    let right = catalog([
        record(20, 0.0, 0.2, 2016.0),
        record(21, 4.0, FRAC_PI_2, 2016.0),
        record(22, 5.0, -FRAC_PI_2, 2016.0),
    ]);
    let matches = left
        .cross_match(&right, options(Radians::new(0.0)))
        .unwrap();
    assert_eq!(matches.len(), 3);
    for (index, candidate) in matches.iter().enumerate() {
        assert_eq!(
            (candidate.left.index, candidate.right.index),
            (index, index)
        );
        assert_eq!(candidate.separation.value(), 0.0);
        assert!(!candidate.left.assumed_stationary && !candidate.right.assumed_stationary);
    }
}

#[test]
fn empty_catalogs_return_no_candidates() {
    let empty = LargeStarCatalog::new();
    let one = catalog([record(1, 1.0, 0.2, 2016.0)]);
    for (left, right) in [(&empty, &one), (&one, &empty), (&empty, &empty)] {
        assert!(left
            .cross_match(right, options(Radians::new(PI)))
            .unwrap()
            .is_empty());
    }
}

#[test]
fn both_catalogs_propagate_from_different_epochs_and_cross_ra_zero() {
    // At the equator pure eastward rectilinear motion changes RA by
    // atan(pm * dt). Opposite offsets at the two epochs meet at RA = 0.001°.
    let rate = Arcseconds::new(8.0).to::<Radian>().value();
    let target_ra = 0.001_f64.to_radians();
    let mut early = record(7, TAU + target_ra - (16.0 * rate).atan(), 0.0, 2000.0);
    let mut late = record(7, target_ra + (10.0 * rate).atan(), 0.0, 2026.0);
    for row in [&mut early, &mut late] {
        row.pm_ra_cosdec = Some(8000.0);
        row.pm_dec = Some(0.0);
    }
    let left = catalog([early]);
    let right = catalog([late]);
    let matches = left
        .cross_match(&right, options(Arcseconds::new(0.01).to()))
        .unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        (matches[0].left.source_id, matches[0].right.source_id),
        (7, 7)
    );
    assert!(matches[0].separation.value() < 2e-15);
}

#[test]
fn eastward_motion_near_and_at_the_pole_is_retained() {
    let rate = Arcseconds::new(10.0).to::<Radian>().value();
    let offset = rate.atan();
    // From the pole, the recorded RA = 0 defines east as the +y direction.
    let right = catalog([record(2, FRAC_PI_2, FRAC_PI_2 - offset, 2016.0)]);
    for dec in [FRAC_PI_2, FRAC_PI_2 - 5e-13] {
        let mut early = record(1, 0.0, dec, 2015.0);
        early.pm_ra_cosdec = Some(10_000.0);
        early.pm_dec = Some(0.0);
        let matches = catalog([early])
            .cross_match(&right, options(Radians::new(1e-12)))
            .unwrap();
        assert_eq!(matches.len(), 1);
        assert!(matches[0].separation.value() < 6e-13);
    }
}

#[test]
fn deterministic_ties_and_ambiguity_counts_include_hidden_pairs() {
    let left = catalog([
        record(5, 0.0, 0.0, 2016.0),
        record(5, 0.0, 0.0, 2016.0),
        record(5, 0.0, 0.2, 2016.0),
    ]);
    // Equal meridional offsets, in the reverse order of the declination index.
    let right = catalog([
        record(5, 0.0, 0.1, 2016.0),
        record(5, 0.0, -0.1, 2016.0),
        record(5, 0.0, 0.2, 2016.0),
    ]);
    let radius = Radians::new(0.11);
    let opts = options(radius);
    assert_eq!(opts.epoch(), JulianYears::new(2016.0));
    assert_eq!(opts.radius(), radius);
    assert_eq!(opts.mode(), CrossMatchMode::AllCandidates);
    assert_eq!(opts.missing_proper_motion(), MissingProperMotion::Error);
    let all = left.cross_match(&right, opts).unwrap();
    let indices: Vec<_> = all.iter().map(|c| (c.left.index, c.right.index)).collect();
    assert_eq!(indices, [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2), (2, 0)]);
    assert_eq!(all[0].separation, all[1].separation);
    assert_eq!(
        (all[0].left_candidate_count, all[0].right_candidate_count),
        (2, 3)
    );
    assert_eq!(left.cross_match(&right, opts).unwrap(), all);
    let opts = opts.with_mode(CrossMatchMode::Nearest);
    assert_eq!(opts.mode(), CrossMatchMode::Nearest);
    let nearest = left.cross_match(&right, opts).unwrap();
    assert_eq!(nearest, [all[0], all[2], all[4]]);
    assert_eq!(nearest[0].right_candidate_count, 3);
}

#[test]
fn missing_motion_policies_distinguish_same_epoch_from_stale_positions() {
    let missing = record(1, 0.0, 0.0, 2000.0);
    let mut partial = missing.clone();
    partial.pm_dec = Some(1e9);
    let left = catalog([record(3, 0.0, 0.0, 2016.0)]);
    for compromised in [missing, partial] {
        let right = catalog([compromised, record(2, 0.0, 0.0, 2016.0)]);
        let opts = options(Radians::zero());
        assert_eq!(
            left.cross_match(&right, opts),
            Err(CrossMatchError::MissingProperMotion {
                side: CrossMatchSide::Right,
                index: 0
            })
        );
        let opts = opts.with_missing_proper_motion(MissingProperMotion::Skip);
        let skipped = left.cross_match(&right, opts).unwrap();
        assert_eq!(skipped.len(), 1);
        assert_eq!(skipped[0].right.index, 1);
        assert!(!skipped[0].right.assumed_stationary);
        let opts = opts.with_missing_proper_motion(MissingProperMotion::AssumeStationary);
        assert_eq!(
            opts.missing_proper_motion(),
            MissingProperMotion::AssumeStationary
        );
        let stationary = left.cross_match(&right, opts).unwrap();
        assert_eq!(stationary.len(), 2);
        assert!(stationary[0].right.assumed_stationary);
        assert!(!stationary[1].right.assumed_stationary);
        assert!(
            right.cross_match(&left, opts).unwrap()[0]
                .left
                .assumed_stationary
        );
    }
}

#[test]
fn invalid_options_and_rows_are_reported_even_without_possible_pairs() {
    let empty = LargeStarCatalog::new();
    for epoch in [f64::NAN, f64::INFINITY] {
        let opts = CrossMatchOptions::new(JulianYears::new(epoch), Radians::zero());
        assert_eq!(
            empty.cross_match(&empty, opts),
            Err(CrossMatchError::InvalidEpoch)
        );
    }
    for radius in [f64::NAN, f64::INFINITY, -0.1, PI + 0.1] {
        assert_eq!(
            empty.cross_match(&empty, options(Radians::new(radius))),
            Err(CrossMatchError::InvalidRadius)
        );
    }

    let row = record(1, 0.0, 0.0, 2000.0);
    let invalid = [
        (
            CatalogRecord {
                ra: Radians::new(-0.1),
                ..row.clone()
            },
            "ra",
        ),
        (
            CatalogRecord {
                ra: Radians::new(TAU + 0.1),
                ..row.clone()
            },
            "ra",
        ),
        (
            CatalogRecord {
                ra: Radians::new(f64::NAN),
                ..row.clone()
            },
            "ra",
        ),
        (
            CatalogRecord {
                dec: Radians::new(FRAC_PI_2 + 0.1),
                ..row.clone()
            },
            "dec",
        ),
        (
            CatalogRecord {
                dec: Radians::new(f64::NEG_INFINITY),
                ..row.clone()
            },
            "dec",
        ),
        (
            CatalogRecord {
                epoch_jyr: f64::NAN,
                ..row.clone()
            },
            "epoch_jyr",
        ),
        (
            CatalogRecord {
                pm_ra_cosdec: Some(f64::NAN),
                ..row.clone()
            },
            "pm_ra_cosdec",
        ),
        (
            CatalogRecord {
                pm_dec: Some(f64::INFINITY),
                ..row.clone()
            },
            "pm_dec",
        ),
        (
            CatalogRecord {
                epoch_jyr: -f64::MAX,
                pm_ra_cosdec: Some(1e20),
                pm_dec: Some(0.0),
                ..row
            },
            "epoch propagation",
        ),
    ];
    let opts = options(Radians::zero()).with_missing_proper_motion(MissingProperMotion::Skip);
    for (bad, field) in invalid {
        assert_eq!(
            catalog([bad]).cross_match(&empty, opts),
            Err(CrossMatchError::InvalidRecord {
                side: CrossMatchSide::Left,
                index: 0,
                field
            })
        );
    }
}

#[test]
fn separation_remains_resolved_from_tiny_angles_to_antipodes() {
    // Near zero the meridional distance is exactly the declination offset.
    let left = catalog([record(0, 0.0, 0.0, 2016.0)]);
    let tiny = 1e-15;
    let right = catalog([
        record(1, 0.0, tiny, 2016.0),
        record(2, PI, 0.0, 2016.0),
        record(3, PI - 1e-9, 0.0, 2016.0),
    ]);
    let matches = left.cross_match(&right, options(Radians::new(PI))).unwrap();
    assert_eq!(matches.len(), 3);
    assert!((matches[0].separation.value() - tiny).abs() < 1e-29);
    assert!((matches[1].separation.value() - (PI - 1e-9)).abs() < 1e-15);
    assert_eq!(matches[2].separation.value(), PI);

    let hits = left
        .cross_match(&right, options(matches[0].separation))
        .unwrap();
    assert_eq!(hits.len(), 1);
    let mut expected = matches[0];
    expected.left_candidate_count = 1;
    assert_eq!(hits[0], expected);
    let hits = left
        .cross_match(&right, options(Radians::new(tiny * 0.9)))
        .unwrap();
    assert!(hits.is_empty());
}

#[test]
fn separations_agree_with_erfa_across_ra_zero_and_near_an_antipode() {
    let left = catalog([record(0, TAU - 5e-6, 1.23456, 2016.0)]);
    let right = catalog([
        record(1, 0.0000123, 1.2345600001, 2016.0),
        record(2, 3.1415873464, -1.23456, 2016.0),
        record(3, 4.321, 0.876, 2016.0),
    ]);
    // ERFA v2.0.1 eraSeps on the input radians, in right row order.
    let reference = [
        5.707_901_762_887_911e-6,
        3.141_592_552_236_662_3,
        0.870_310_991_890_170_4,
    ];
    let hits = left.cross_match(&right, options(Radians::new(PI))).unwrap();
    assert_eq!(hits.len(), reference.len());
    for hit in hits {
        assert!((hit.separation.value() - reference[hit.right.index]).abs() < 2e-15);
    }
}

#[test]
fn both_motion_components_follow_the_expected_initial_bearing() {
    let ra = 1.0;
    let dec: f64 = 1.2;
    let mut early = record(1, ra, dec, 1966.0);
    early.pm_ra_cosdec = Some(6000.0);
    early.pm_dec = Some(-8000.0);

    // Solve the spherical destination problem using a 3:4:5 motion triangle.
    // The direction of the projected point has traveled atan(speed * dt)
    // along a great circle with a bearing measured eastward from north.
    let (sin_arc, cos_arc) = Arcseconds::new(500.0)
        .to::<Radian>()
        .value()
        .atan()
        .sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();
    let target_dec = (sin_dec * cos_arc - cos_dec * sin_arc * 0.8).asin();
    let delta_ra = (0.6 * sin_arc * cos_dec).atan2(cos_arc - sin_dec * target_dec.sin());
    let right = catalog([record(2, ra + delta_ra, target_dec, 2016.0)]);
    let hits = catalog([early])
        .cross_match(&right, options(Radians::new(1e-12)))
        .unwrap();
    assert_eq!(hits.len(), 1);
    assert!(hits[0].separation.value() < 2e-15);
}
