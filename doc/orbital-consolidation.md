# Orbital API consolidation and coordinated release

Resolves Siderust #113 with upstream affn #36, keplerian #9 and tempoch #23.
Siderust remains **0.15.0**. Its PRs #114 and #112 are not yet merged; the
upstream patch releases are available from crates.io.

## Baseline and ownership inventory

The baseline was fetched before modification; all target working trees were clean.
PR #112 was fast-forwarded from b776b8a to its latest remote 0464b2f, preserving
seven existing commits. Its baseline CI passed all 11 jobs. The integration branch
starts from latest main **57554f2**, not from the SGP4 branch.

| Crate | Baseline manifest | Baseline commit | Integrated version |
|---|---|---|---|
| affn | 0.12.0 | 8c60a50 | 0.12.1 |
| affn-derive | 0.12.0 | same repository | unchanged |
| keplerian | 0.5.0 | e751482 | 0.5.1 |
| tempoch / tempoch-core / tempoch-ffi | 0.11.0 | 12d96e8 | 0.11.1 |
| tempoch-validation | 0.1.0 | same repository | unchanged (not published) |
| siderust / siderust-ffi | 0.15.0 | 57554f2 / PR112 0464b2f | 0.15.0 |
| qtty | 0.8.6 | 761cfc0 | unchanged |
| principia | 0.7.0 | 7527ccc | unchanged |

Tags/releases, open PRs, full issue descriptions/discussions, READMEs, manifests,
changelogs, feature matrices, publication workflows and CI/coverage gates were
inspected. There were no applicable AGENTS.md or separate contribution files in
the target repositories. Existing release tags were not recreated. Compatible
patch releases preserve the current 0.x dependency type identities across
principia, optica, affn and the time crates.

| Existing implementation | Owner / consumers | Decision, reduction and risks |
|---|---|---|
| Siderust `astro/conic_equations.rs` elliptic/hyperbolic true-anomaly formulas; keplerian `anomaly.rs` | keplerian; Siderust orbit/conic propagation | Real algorithm duplication: delegate both to canonical typed functions. Elliptic canonical output is [0,2π); positions are equivalent to the old signed branch. Test wraps, negative/large angles, circular and near-parabolic limits. |
| keplerian `elements.rs::rotate_pqw`; Siderust `conic_equations.rs` rotations and `orbit.rs::OrientationTrig` | affn `Rotation3`; both consumers | Remove independently maintained orbital matrix/trig implementations. Cache one affn matrix in PreparedOrbit, compute it once per keplerian Cartesian conversion. Existing numerical/FFI oracles remain unchanged. |
| keplerian public signed-axis/radian elements vs affn `OrientedConic`, `ConicOrientation`, `SemiMajorAxisParam` | affn geometry; keplerian dynamics | Preserve public mutable fields and Serde; add a checked, on-demand projection. Positive hyperbolic magnitude, degree units and exact parabolic classification differ from legacy Keplerian conventions; never silently substitute validation/tolerance policies. No duplicate stored representation. |
| affn Cartesian Position/Velocity, keplerian CartesianState, principia DynamicsState, Siderust TemeState | affn geometry, domain aggregates in their owners | Storage similarity, no independent vector algorithm. Keep existing aggregates and public field/epoch contracts. Siderust OrbitState already aliases principia with TT; UTC SGP4 cannot use a continuous-scale state. No extra wrapper. |
| affn `BlockDiagRotation6`, Siderust TEME rotation providers | affn instantaneous geometry; Siderust astronomical conventions | Correct documentation, retain diag(R,R). A coordinate-state/covariance change needs Rdot and possibly origin velocity. No actual full-state caller supplies that derivative, so no unused kinematic framework or untested TEME state conversion is introduced. |
| SGP4 `time.rs::minutes_between`, `offset_epoch`; Tempoch calendar bridge / ExactDuration | Tempoch time conversion, Siderust model adapter | Existing precise explicit-context APIs and chrono already own the algorithms. Document/test reuse; retain short model/error adapters and 1,440 nominal minutes per UTC civil day. Do not use SI subtraction or subtract rounded JD scalars. |
| Elliptic radius `a*(1-e*cos(E))` | Single Siderust wrapper | One expression with one independent downstream consumer does not justify a new public helper. It remains local, unit convention unchanged. |
| SGP4 WGS constants, AFSPC/IAU sidereal initialization, un-Kozai/resonance/drag kernel | SGP4 backend; Siderust orchestration pending corrected backend | Scientifically distinct from two-body dynamics and modern IAU frame code. Retain unchanged; never transfer to keplerian/principia or replace with generic GM constants. |

There are no new dependency cycles. QTTY continues to own dimensions, constants
and units; affn has no epochs/EOP; keplerian has no astronomical time/TLE; Tempoch
has no SGP4 model; principia retains continuous-time integration. Public orbital
names, argument/return types, coordinate wrappers, Serde formats and FFI layouts
remain compatible. `PreparedOrbit`'s private cache grows from six to nine f64s
(+24 bytes) in exchange for removing the independent rotation implementation;
it is not an FFI layout. Tempoch's version-reporting C symbol reports 1101 instead
of 1100, with no symbol/layout changes. Tests explicitly check the new version.

Tempoch full parallel validation exposed an existing test-only synthetic-table
race. Thread-local overrides and unwind restoration isolate the fixtures without
changing production runtime-data behavior or serializing/disabling the suite.

## Registry dependencies and integration order

The upstream PRs affn #37, keplerian #10 and Tempoch #24 have merged.
The integration uses published affn **0.12.1**, keplerian **0.5.1**, and
tempoch / tempoch-core / tempoch-ffi **0.11.1**. Cargo resolves all of these
directly from crates.io; `Cargo.lock` records their registry checksums.
There is no root `[patch.crates-io]` override, local dependency checkout, or
CI-only dependency substitution. Compatible transitive consumers continue
to unify on one version of each crate family.

1. Merge Siderust PR #114 before PR #112. The SGP4 branch already contains
   this refactor in its ancestry, so refresh PR #112 against main after #114
   is squash-merged and confirm its effective diff.
2. Verify registry resolution and rerun the complete Rust/FFI, numerical,
   no_std/Serde/embedded and coverage suites before releasing Siderust 0.15.0.
   The native-versus-official SGP4 backend decision remains separate in #98.

To inspect the resolved graph and run the main local checks:

```sh
cargo tree -d
cargo tree --workspace -i affn
cargo tree --workspace -i keplerian
cargo tree --workspace -i tempoch-core
cargo tree --workspace -i tempoch-ffi
cargo metadata --format-version 1 --locked > registry-metadata.json
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo check --workspace --all-features
SIDERUST_JPL_STUB=all cargo test --workspace --all-features
SIDERUST_JPL_STUB=all cargo test --doc --workspace --all-features
```

The metadata must show registry sources and the versions above, without
duplicate versions of these families. The JPL stub matches normal CI and
does not alter scientific reference fixtures.

## Maintained LOC

Measure with `uv run scripts/measure-orbital-loc.py /path/to/repo REF`. The script
pins its Rust syntax parser, uses only tracked source at the selected revision,
excludes cfg(test) items/statements and comments, and counts nonblank physical
Rust code lines in package src directories. Tests, text fixtures, generated files,
headers and Markdown are excluded and reported separately. Build scripts,
benchmarks, examples and measurement/release scripts are not production src.
This is a consistent physical-line measure, not a compression or coverage metric.

| Repository | Production before | Production after | Net |
|---|---:|---:|---:|
| affn | 7,166 | 7,166 | 0 |
| keplerian | 1,629 | 1,624 | -5 |
| tempoch | 5,662 | 5,665 | +3 |
| siderust (#113 vs main) | 39,433 | 39,395 | -38 |
| **Total** | **53,890** | **53,850** | **-40** |

Keplerian's checked geometry adapter adds code while deleting its rotation
implementation; Tempoch's test-isolation cfg branches add three physical
production lines without a new production algorithm. Siderust removes its
independent orientation cache and two anomaly formulas. Affn adds no production
code. This modest net reduction is secondary to eliminating independently
maintained mathematical implementations. No validation or reference data was
removed. The single radius expression, domain aggregates and astronomy/SGP4
specializations remain for the concrete reasons in the inventory above.

Excluded raw-line counts (comments inside production Rust are additionally
excluded; Markdown figures include changelogs/design documents):

| Repository | Test lines before → after | Fixture lines before → after | Generated/header lines before → after | Markdown lines before → after |
|---|---:|---:|---:|---:|
| affn | 6,895 → 6,952 | 0 → 0 | 0 → 0 | 1,288 → 1,338 |
| keplerian | 1,522 → 1,620 | 0 → 0 | 0 → 0 | 169 → 207 |
| tempoch | 7,480 → 7,576 | 0 → 0 | 612 → 612 | 1,110 → 1,168 |
| siderust | 36,831 → 36,888 | 1,695 → 1,695 | 2,006 → 2,006 | 3,697 → 3,921 |

## Validation

Upstream local validation includes workspace all-target/all-feature Clippy,
check/tests/doctests, formatting, no_std/alloc/Serde and supported thumb builds.
Unchanged local line-coverage gates pass: affn **91.01%** (85% required),
keplerian **93.60%** (90%), Tempoch **94.41%** (90%). GitHub Actions additionally
retains keplerian audit/deny and Miri; no test or threshold is disabled.

Siderust validates the registry-resolved dependency graph with all-feature workspace
checks, tests and doctests, all-target Clippy, no_std/Serde/embedded builds and
unchanged CI feature/coverage jobs. New independent comparisons preserve legacy
anomaly positions and rotations at wraps, circular/high-e/near-parabolic limits,
negative/large angles and scaled lengths. Existing tests cover large propagation
intervals, hyperbolic orbits, catalogs, frames/centers, Serde and the C ABI.
All local runtime tests, 310 doctests (26 pre-existing ignores), and the
no_std/Serde/embedded matrix pass. A transient /tmp-backed linker failure was
resolved by using disk-backed TMPDIR and two doctest threads; no source/test
change was needed. CI check links and final results are in the PR descriptions.

## Remaining risks and follow-ups

The official sgp4 2.4.0 G520 boundary discrepancy is a separate #98 backend
migration concern. The validated native kernel stays in PR #112; no generic
solver replaces it and none of its mathematics moves to another ecosystem crate.
A real rotating-frame state conversion still needs tested astronomical Rdot and
origin-velocity data; the documentation now states this limitation explicitly.
Neither item blocks reviewing this generic API consolidation. Merge #114 before #112, then complete final regression checks before release.
