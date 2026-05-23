// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Synthetic POD smoke test.
//!
//! Exercises the full domain-driven public surface of `siderust-pod` —
//! `problem`, `run`, `estimation` — using only in-memory synthetic data,
//! without any external files, TLE strings, or real ephemeris datasets.
//!
//! The test is intentionally small: it verifies that the typed module
//! boundaries are coherent and that a minimal WLS solve round-trips
//! without panicking or returning an unexpected error variant.

use qtty::Second;
use tempoch::{J2000Seconds, TT};

use siderust::pod::estimation::{NormalEquations, WlsSolverError};
use siderust::pod::problem::arc::{ArcDefinition, ArcId};
use siderust::pod::problem::covariance::ParameterCovariance;
use siderust::pod::problem::parameter::{Parameter, ParameterKind, ParameterOrdering};
use siderust::pod::run::dataset::DatasetRef;
use siderust::pod::run::manifest::RunManifest;

fn tt(secs: f64) -> tempoch::Time<TT> {
    J2000Seconds::<TT>::try_new(Second::new(secs)).unwrap()
}

// ── Problem layer ─────────────────────────────────────────────────────────────

#[test]
fn problem_arc_roundtrips_display() {
    let id = ArcId::new("arc-001");
    assert_eq!(id.as_str(), "arc-001");
    assert!(format!("{id:?}").contains("arc-001"));
}

#[test]
fn problem_arc_definition_stores_fields() {
    let id = ArcId::new("arc-A");
    let def = ArcDefinition {
        id: id.clone(),
        start: tt(0.0),
        stop: tt(86_400.0),
        step_hint: Some(Second::new(60.0)),
    };
    assert_eq!(def.id, id);
    assert!(def.is_valid());
    assert!((def.duration().value() - 86_400.0).abs() < 1e-6);
}

#[test]
fn problem_parameter_ordering_maps_kind_to_index() {
    let kinds = [
        ParameterKind::StatePositionX,
        ParameterKind::StatePositionY,
        ParameterKind::StatePositionZ,
        ParameterKind::StateVelocityX,
        ParameterKind::StateVelocityY,
        ParameterKind::StateVelocityZ,
    ];
    let ordering = ParameterOrdering {
        params: kinds
            .iter()
            .cloned()
            .map(|k| Parameter {
                kind: k,
                initial_value: 0.0,
                apriori_sigma: None,
            })
            .collect(),
    };
    for (i, k) in kinds.iter().enumerate() {
        assert_eq!(ordering.index_of(k), Some(i));
    }
    assert_eq!(ordering.len(), 6);
}

#[test]
fn problem_parameter_round_trips_apriori() {
    let p = Parameter {
        kind: ParameterKind::DragScale,
        initial_value: 1.0,
        apriori_sigma: Some(0.1),
    };
    assert_eq!(p.kind, ParameterKind::DragScale);
    assert!((p.initial_value - 1.0).abs() < f64::EPSILON);
    assert!((p.apriori_sigma.unwrap() - 0.1).abs() < f64::EPSILON);
}

#[test]
fn problem_covariance_diagonal_roundtrip() {
    let n = 3usize;
    let sigmas = [1.0_f64, 2.0, 3.0];
    let mut data = vec![0.0f64; n * n];
    for (i, &s) in sigmas.iter().enumerate() {
        data[i * n + i] = s * s;
    }
    let ordering = ParameterOrdering {
        params: sigmas
            .iter()
            .map(|&s| Parameter {
                kind: ParameterKind::DragScale,
                initial_value: 0.0,
                apriori_sigma: Some(s),
            })
            .collect(),
    };
    let cov = ParameterCovariance {
        params: ordering,
        data,
    };
    assert_eq!(cov.n(), n);
    for (i, &s) in sigmas.iter().enumerate() {
        let diag = cov.data[i * n + i];
        assert!(
            (diag - s * s).abs() < 1e-12,
            "variance[{i}] = {diag}, expected {}",
            s * s
        );
    }
    cov.validate_symmetric(1e-12).unwrap();
    cov.validate_diagonal_nonneg().unwrap();
}

// ── Run layer ────────────────────────────────────────────────────────────────

#[test]
fn run_dataset_ref_stores_kind_and_path() {
    let d = DatasetRef::from_bytes("/data/igs/igs12345.sp3", "sp3", b"dummy");
    assert_eq!(d.kind, "sp3");
    assert!(d.path.to_str().unwrap().ends_with(".sp3"));
    assert!(!d.sha256.is_empty());
    assert_eq!(d.bytes, 5);
}

#[test]
fn run_manifest_canonicalize_sorts_io() {
    let mut manifest = RunManifest {
        run_id: "run-001".to_string(),
        tool_version: "siderust-pod 0.0.0".to_string(),
        config_sha256: "deadbeef".to_string(),
        inputs: vec![
            DatasetRef::from_bytes("b.sp3", "sp3", b"x"),
            DatasetRef::from_bytes("a.rnx", "rinex", b"x"),
        ],
        outputs: vec![
            DatasetRef::from_bytes("z.csv", "residuals", b"x"),
            DatasetRef::from_bytes("a.sp3", "orbit", b"x"),
        ],
        started_at: "2024-01-01T00:00:00Z".to_string(),
        finished_at: "2024-01-01T01:00:00Z".to_string(),
    };
    manifest.canonicalize();
    // Inputs: rinex < sp3 alphabetically on kind
    assert_eq!(manifest.inputs[0].kind, "rinex");
    assert_eq!(manifest.inputs[1].kind, "sp3");
    // Outputs: orbit < residuals alphabetically on kind
    assert_eq!(manifest.outputs[0].kind, "orbit");
    assert_eq!(manifest.outputs[1].kind, "residuals");
}

// ── Estimation layer (WLS) ───────────────────────────────────────────────────

/// Build a trivial 2-parameter normal equations system.
///
/// True state: [x=1.0, y=2.0].
/// Observations: `x = 1.0 ± 0.1`, `y = 2.0 ± 0.2`, `x+y = 3.0 ± 0.15`.
fn build_synthetic_normal_equations() -> NormalEquations {
    let mut neq = NormalEquations::new(2);
    neq.add_row(&[(0, 1.0)], 1.0, 0.1).unwrap();
    neq.add_row(&[(1, 1.0)], 2.0, 0.2).unwrap();
    neq.add_row(&[(0, 1.0), (1, 1.0)], 3.0, 0.15).unwrap();
    neq
}

#[test]
fn wls_solve_overdetermined_system_converges() {
    let neq = build_synthetic_normal_equations();
    let result = neq
        .solve()
        .expect("WLS solve should succeed on well-conditioned system");
    assert!(
        (result.update[0] - 1.0).abs() < 1e-6,
        "x component: expected ~1.0, got {}",
        result.update[0]
    );
    assert!(
        (result.update[1] - 2.0).abs() < 1e-6,
        "y component: expected ~2.0, got {}",
        result.update[1]
    );
}

#[test]
fn wls_reduced_chi2_is_finite_and_positive() {
    let neq = build_synthetic_normal_equations();
    let result = neq.solve().unwrap();
    let chi2 = result.reduced_chi2();
    assert!(chi2.is_finite(), "reduced chi2 should be finite");
    assert!(chi2 >= 0.0, "reduced chi2 should be non-negative");
}

#[test]
fn wls_underdetermined_system_errors_gracefully() {
    // Zero observations → singular normal matrix → should return an error
    let neq = NormalEquations::new(2);
    let result = neq.solve();
    assert!(
        matches!(result, Err(WlsSolverError::NotPositiveDefinite(_))),
        "Expected NotPositiveDefinite error, got: {result:?}"
    );
}

// ── End-to-end: problem → run → estimation ───────────────────────────────────

/// Synthetic single-arc POD cycle:
/// 1. Define the arc and parameter ordering.
/// 2. Record run metadata.
/// 3. Build normal equations from simulated residuals.
/// 4. Solve and assert convergence.
#[test]
fn synthetic_pod_cycle_problem_run_estimation() {
    // 1. Problem
    let arc = ArcDefinition {
        id: ArcId::new("synthetic-arc-001"),
        start: tt(0.0),
        stop: tt(3600.0),
        step_hint: Some(Second::new(30.0)),
    };
    assert!(arc.is_valid());

    let ordering = ParameterOrdering {
        params: vec![
            Parameter {
                kind: ParameterKind::StatePositionX,
                initial_value: 6378.0,
                apriori_sigma: Some(1.0),
            },
            Parameter {
                kind: ParameterKind::StatePositionY,
                initial_value: 0.0,
                apriori_sigma: Some(1.0),
            },
            Parameter {
                kind: ParameterKind::StatePositionZ,
                initial_value: 0.0,
                apriori_sigma: Some(1.0),
            },
        ],
    };
    assert_eq!(ordering.len(), 3);
    assert_eq!(ordering.index_of(&ParameterKind::StatePositionZ), Some(2));

    // 2. Run manifest stub
    let manifest = RunManifest {
        run_id: format!("{}-run", arc.id.as_str()),
        tool_version: "siderust-pod 0.0.0".to_string(),
        config_sha256: "0000".to_string(),
        inputs: vec![],
        outputs: vec![],
        started_at: "2024-01-01T00:00:00Z".to_string(),
        finished_at: "2024-01-01T00:00:01Z".to_string(),
    };
    assert!(manifest.run_id.contains("synthetic"));

    // 3. Simulated measurements: position corrections [+0.5, -0.3, +0.1] km
    let mut neq = NormalEquations::new(3);
    let sigma = 0.01_f64; // 10 m position noise
    neq.add_row(&[(0, 1.0)], 0.5, sigma).unwrap();
    neq.add_row(&[(1, 1.0)], -0.3, sigma).unwrap();
    neq.add_row(&[(2, 1.0)], 0.1, sigma).unwrap();
    // Extra redundant observations for overdetermination
    neq.add_row(&[(0, 1.0), (1, 1.0)], 0.2, sigma * 2.0_f64.sqrt())
        .unwrap();
    neq.add_row(&[(1, 1.0), (2, 1.0)], -0.2, sigma * 2.0_f64.sqrt())
        .unwrap();

    // 4. Solve
    let result = neq.solve().expect("synthetic POD solve must succeed");
    assert!(
        (result.update[0] - 0.5).abs() < 0.01,
        "Δx ≈ 0.5 km, got {}",
        result.update[0]
    );
    assert!(
        (result.update[1] - (-0.3)).abs() < 0.01,
        "Δy ≈ -0.3 km, got {}",
        result.update[1]
    );
    assert!(
        (result.update[2] - 0.1).abs() < 0.01,
        "Δz ≈ 0.1 km, got {}",
        result.update[2]
    );
    assert!(result.reduced_chi2().is_finite());
}
