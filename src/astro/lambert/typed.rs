//! Typed entry-points for Lambert's problem.
//!
//! These functions are the **public API** of the crate. They operate on
//! `affn` typed positions and `affn` typed velocities, with the time of
//! flight expressed as a `qtty::Second` (chosen over
//! [`tempoch::Period`](https://docs.rs/tempoch) because Lambert ToF is
//! semantically a *duration*, not a half-open instant interval — see the
//! design note at the bottom of this module).
//!
//! For interoperability with raw-array callers (mission-design search
//! loops, FFI bridges, ...) the numeric kernel is exposed via
//! [`crate::solve_lambert`] / [`crate::solve_lambert_n_rev`].

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use qtty::length::Kilometer;
use qtty::Second;

use super::error::LambertError;
use super::izzo::{
    solve_lambert as solve_lambert_arr, solve_lambert_n_rev as solve_lambert_n_rev_arr,
    LambertBranch, LambertDiagnostics, LambertSolution, NRevBranch,
};

/// Typed Lambert solution — departure / arrival velocities plus diagnostics.
///
/// Velocities are tagged with the same reference frame `F` as the input
/// positions. They are *free vectors*, so no [`ReferenceCenter`] tag is
/// attached (cf. `affn::cartesian::Velocity = Vector<F, U>`).
///
/// # Examples
///
/// ```
/// use siderust_pod::lambert::TypedLambertSolution;
/// // The struct just bundles the two velocities and diagnostics; see
/// // [`siderust_pod::lambert::lambert`] for end-to-end usage.
/// fn assert_solution_layout<F: affn::frames::ReferenceFrame>(
///     sol: TypedLambertSolution<F>,
/// ) {
///     let _v1 = sol.v1;
///     let _v2 = sol.v2;
///     let _diag = sol.diagnostics;
/// }
/// ```
#[derive(Debug, Clone, Copy)]
pub struct TypedLambertSolution<F: ReferenceFrame> {
    /// Departure velocity at `r1`, in the same frame as the input positions.
    pub v1: Velocity<F, KmPerSecond>,
    /// Arrival velocity at `r2`, in the same frame as the input positions.
    pub v2: Velocity<F, KmPerSecond>,
    /// Householder iteration diagnostics from the underlying numeric solver.
    pub diagnostics: LambertDiagnostics,
}

/// Solve Lambert's problem (single revolution) on typed inputs.
///
/// `r1`, `r2` are typed [`Position`]s in the same reference frame `F`
/// and reference center `C`; `tof` is a [`qtty::Second`] duration; `mu`
/// is a typed [`GravitationalParameter`] (km³/s²). The function returns
/// the corresponding free-vector velocities at `r1` and `r2`.
///
/// The two endpoints **must** share both the frame and the center —
/// this is enforced statically through the type signature, mirroring
/// the affine-space invariant that "subtracting two `Position`s only
/// makes sense within a single frame and centre".
///
/// # Errors
///
/// Wraps the underlying [`LambertError`] taxonomy unchanged.
///
/// # Examples
///
/// ```
/// // no extra center import — () implements ReferenceCenter
/// use affn::frames::ICRS;
/// use affn::cartesian::Position;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometer;
/// use qtty::{Quantity, Second};
/// use siderust_pod::lambert::{lambert, LambertBranch};
///
/// let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
/// let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
/// let tof = Second::new(4_560.0);
/// let mu = GravitationalParameter::new(398_600.4418);
///
/// let sol = lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
/// assert!((sol.v1.x().value() - 2.058913).abs() < 1e-3);
/// ```
pub fn lambert<C, F>(
    r1: Position<C, F, Kilometer>,
    r2: Position<C, F, Kilometer>,
    tof: Second,
    mu: GravitationalParameter,
    branch: LambertBranch,
) -> Result<TypedLambertSolution<F>, LambertError>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let LambertSolution {
        v1,
        v2,
        diagnostics,
    } = solve_lambert_arr(
        position_to_array(&r1),
        position_to_array(&r2),
        tof.value(),
        mu.value(),
        branch,
    )?;

    Ok(TypedLambertSolution {
        v1: Velocity::<F, KmPerSecond>::new(v1[0], v1[1], v1[2]),
        v2: Velocity::<F, KmPerSecond>::new(v2[0], v2[1], v2[2]),
        diagnostics,
    })
}

/// Solve Lambert's problem with `N ≥ 1` complete revolutions, on typed inputs.
///
/// `side` selects between the *low-energy* (`Left`) and *high-energy*
/// (`Right`) Householder roots that bracket the time-of-flight minimum
/// for this `N`.
///
/// # Errors
///
/// Wraps the underlying [`LambertError`] variants. In particular,
/// requesting `N` beyond the geometrically admissible maximum returns
/// [`LambertError::RevolutionsExceedNMax`] with the maximum that *is*
/// admissible at the supplied `tof`.
///
/// # Examples
///
/// ```
/// // no extra center import — () implements ReferenceCenter
/// use affn::frames::ICRS;
/// use affn::cartesian::Position;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometer;
/// use qtty::Second;
/// use siderust_pod::lambert::{lambert_n_rev, LambertBranch, NRevBranch};
///
/// let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
/// let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
/// let tof = Second::new(21_600.0); // 6 h
/// let mu = GravitationalParameter::new(398_600.4418);
///
/// let sol = lambert_n_rev(
///     r1, r2, tof, mu,
///     LambertBranch::Prograde, 1, NRevBranch::Left,
/// )
/// .unwrap();
/// assert_eq!(sol.diagnostics.revolutions, 1);
/// ```
pub fn lambert_n_rev<C, F>(
    r1: Position<C, F, Kilometer>,
    r2: Position<C, F, Kilometer>,
    tof: Second,
    mu: GravitationalParameter,
    branch: LambertBranch,
    revolutions: u32,
    side: NRevBranch,
) -> Result<TypedLambertSolution<F>, LambertError>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let LambertSolution {
        v1,
        v2,
        diagnostics,
    } = solve_lambert_n_rev_arr(
        position_to_array(&r1),
        position_to_array(&r2),
        tof.value(),
        mu.value(),
        branch,
        revolutions,
        side,
    )?;

    Ok(TypedLambertSolution {
        v1: Velocity::<F, KmPerSecond>::new(v1[0], v1[1], v1[2]),
        v2: Velocity::<F, KmPerSecond>::new(v2[0], v2[1], v2[2]),
        diagnostics,
    })
}

fn position_to_array<C, F>(p: &Position<C, F, Kilometer>) -> [f64; 3]
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    [p.x().value(), p.y().value(), p.z().value()]
}

// ─────────────────────────────────────────────────────────────────────────────
// Design note: Period vs. Second for the time-of-flight argument
// ─────────────────────────────────────────────────────────────────────────────
//
// `tempoch::Period<S> = Interval<Time<S>>` is a *half-open instant
// interval*, i.e. `[start, end)` on a particular time scale. Lambert's
// problem is invariant under translation of the boundary epochs: only
// the duration `Δt = end − start` enters the dynamics. Forcing callers
// to attach a time scale to that duration would (a) require a scale
// choice the algorithm does not need, and (b) break callers that work
// in elapsed-time form (mission-design grids, porkchop searches).
//
// Conversely, `qtty::Second` is exactly the dimensioned duration the
// algorithm needs: a `Quantity<unit::Second>` with no axis tag. Callers
// holding a `Period<S>` can trivially convert via `period.end -
// period.start`, which already yields a `Second` per
// `tempoch::time::Time<S>::Sub`.

#[cfg(test)]
mod tests {
    use super::*;
    // no extra center import — () implements ReferenceCenter
    use affn::frames::ICRS;

    #[test]
    fn typed_zero_rev_matches_array_kernel() {
        let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
        let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
        let tof = Second::new(4_560.0);
        let mu = GravitationalParameter::new(398_600.441_8);

        let typed = lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
        let raw = solve_lambert_arr(
            [15945.34, 0.0, 0.0],
            [12214.83899, 10249.46731, 0.0],
            4_560.0,
            398_600.441_8,
            LambertBranch::Prograde,
        )
        .unwrap();
        for i in 0..3 {
            assert!((typed.v1.as_array()[i].value() - raw.v1[i]).abs() < 1e-12);
            assert!((typed.v2.as_array()[i].value() - raw.v2[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn typed_n_rev_propagates_error() {
        let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
        let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
        let tof = Second::new(4_560.0);
        let mu = GravitationalParameter::new(398_600.441_8);
        let err = lambert_n_rev(
            r1,
            r2,
            tof,
            mu,
            LambertBranch::Prograde,
            1,
            NRevBranch::Left,
        )
        .unwrap_err();
        match err {
            LambertError::RevolutionsExceedNMax { requested, .. } => {
                assert_eq!(requested, 1);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
