# Upstream SGP4 migration audit — blocked

Audited 2026-10-19 against Siderust PR #112, baseline `ab9b7fe135cec4c75e7bb7ebba511c73b3a30343`,
and the latest stable registry `sgp4` **2.4.0**. **The production migration is not
complete. The validated custom core remains in place.** No reference fixtures,
tolerances, public API, FFI, release versions, or CI gates were changed.

## Blocking numerical difference

Upstream's half-day resonance initializer selects the wrong `G520` polynomial
at **exactly eccentricity 0.715**. Its `src/deep_space.rs:496` uses
`eccentricity < 0.715`; the Vallado C++ reference uses the lower polynomial
through equality (`<= 0.715`). This is a numerical backend behavior, not a
unit, epoch, gravity-constant, or frame-convention difference.

Sources:

- [Selected upstream implementation](https://github.com/neuromorphicsystems/sgp4/blob/66b63188adf59d3351abd3d081db4683270f2d54/src/deep_space.rs#L489-L503).
- [Independent C++ reference](https://github.com/brandon-rhodes/python-sgp4/blob/2.25/extension/SGP4.cpp) (`dsinit`, `g520`, `em > 0.715`).
- Unchanged `tests/fixtures/sgp4/synthetic.txt`, generated before this audit;
  provenance and reproduction are in the fixture README.

Five of 105 synthetic states fail in each gravity mode (15 / 315 total).
All failures have eccentricity exactly 0.715 and mean motion 2 rev/day;
10080-minute requests expose the largest differences. Limits remain strictly
**<1e-6 km** and **<1e-9 km/s**, component by component.

| Mode | Failed synthetic states | Max position error (km) | Max velocity error (km/s) |
|---|---:|---:|---:|
| WGS-72 / AFSPC | 5 / 105 | 0.110732743881 | 2.38193318965e-5 |
| WGS-72 / IAU | 5 / 105 | 0.110732762447 | 2.38193358904e-5 |
| WGS-84 / IAU | 5 / 105 | 0.110733047583 | 2.38194970212e-5 |

The complete 33-satellite corpus (666 successful states, two expected errors)
and the 48 independent IAU vectors pass through the experimental typed adapter.
The remaining 300 synthetic vectors pass. This is insufficient for acceptance.
The failure also reproduces directly against upstream with both `std` and
`libm`, independently of Siderust, QTTY, Affn, or Tempoch.

Upstream stores resonance coefficients privately. Gravity selection, sidereal
functions, AFSPC propagation, and fresh resonance state cannot change this
polynomial branch. Changing the input eccentricity would change the user's
scientific input; it is not an acceptable configuration fix.

## Reproduce without changing production dependencies

The standalone diagnostic has its own Cargo workspace and pinned lockfile:

```sh
cargo run --locked --manifest-path scripts/sgp4-backend-audit/Cargo.toml
cargo run --locked --manifest-path scripts/sgp4-backend-audit/Cargo.toml \
  --no-default-features --features libm
```

Both commands deliberately exit **1**, printing all 15 discrepancies and the
per-model maxima. They read the original independent fixtures directly. They
are diagnostic commands, not replacement validation tests or CI exclusions.

As a causal experiment, an isolated copy of the downloaded backend under
`/tmp` was changed only from `< 0.715` to `<= 0.715`. With that diagnostic
patch, all 315 synthetic vectors pass. Maximum component errors become:

| Mode | Position (km) | Velocity (km/s) |
|---|---:|---:|
| WGS-72 / AFSPC | 3.33894e-8 | 8.92086e-12 |
| WGS-72 / IAU | 6.54291e-9 | 1.74616e-12 |
| WGS-84 / IAU | 5.43333e-9 | 1.45395e-12 |

With the adapter's legacy extreme-input range guard retained, all nine typed
numerical tests and both pipeline tests pass against that isolated correction.
**No corrected fork, vendor copy, registry-cache edit, Cargo patch, or changed
orbital input is shipped.** A corrected upstream release is the prerequisite
for completing the requested migration. This audit does not open an upstream
issue or PR.

## Adapter design investigated

- Use `Orbit::from_kozai_elements` with explicit WGS-72 or WGS-84, followed by
  `Constants::new` with the matching sidereal function. Use
  `propagate_afspc_compatibility_mode` only for WGS-72/AFSPC; use `propagate`
  for both IAU modes. Never use the convenience WGS-84 initializer for WGS-72/IAU.
- Keep backend types private, immutable, cloneable, `Send + Sync`. Successful
  propagation needs no heap allocation and fresh local resonance integration
  preserves request-order independence. Avoid the stateful advanced API.
- Use `default-features = false, features = ["libm"]` on every platform.
  Upstream makes `std` and `libm` mutually exclusive, so additive Cargo features
  cannot forward Siderust's `std` feature while keeping unconditional `libm`.
  Upstream allocation and serde features are unnecessary. A downstream graph
  enabling upstream `sgp4/std` would need coordination with this choice.
- Convert mean motion through QTTY `AngularRate<Radian, Minute>`; the unchanged
  independent vectors confirm rounding is within tolerance outside the backend
  boundary bug. Preserve the old finite-input range check: QTTY conversion can
  avoid overflow for `f64::MAX` rev/day that the previous initializer rejected.
  Upstream positivity checks alone also do not reject all NaNs after overflow.
- Reuse `qtty::KmPerSecond` behind the existing public velocity unit alias.
  Keep Affn's existing TEME position and velocity types and transformations.
- Preserve Tempoch validation, nominal UTC calendar differences, leap-label
  rejection, nanosecond range checks, and exact returned target epochs.
  Keep `Minutes` typed until the backend call. Use the existing split calendar
  initialization epoch, retaining AFSPC's legacy scalar-JD rounding, then
  upstream sidereal functions. The backend calendar helpers use a restricted
  calendar formula, including unsigned year arithmetic in AFSPC mode.
- Clean constructors can accept an owned `TimeContext` with or without an
  explicit model. Existing constructors create and own their snapshot.
- Map upstream initialization failures to `InvalidElements`, divergence to
  `Propagation`, and preserve time/epoch variants. Keep non-finite Cartesian
  output checks and diagnostic variant names. FFI remains a Siderust boundary.

An experimental adapter patch and raw local logs are saved in the ignored
`artifacts/sgp4-backend-audit/` directory for resuming the work. The patch is
an investigation snapshot, not a completed or accepted migration.

## Performance experiment

Criterion public-API benchmarks are retained in `benches/sgp4.rs`. Run on an
Intel Core Ultra 9 185H, x86_64 Linux, optimized stable Rust, 30 samples,
1-second warm-up and 2-second measurement windows. Both implementations use
identical benchmark code and WGS-72/AFSPC. Native baseline preceded replacement;
upstream measurements used the unmodified registry backend with `libm`.

```sh
cargo bench --bench sgp4 -- --save-baseline native
# After implementing an acceptable upstream adapter:
cargo bench --bench sgp4 -- --baseline native
```

Values below are Criterion central estimates. Repeated requests use offsets
`i * 60.25` minutes for `i = 0..100`; single requests use 1440.25 minutes.

| Public API operation | Native | Experimental upstream |
|---|---:|---:|
| LEO 06251 initialization | 640.92 ns | 618.60 ns |
| LEO single propagation | 511.76 ns | 513.47 ns |
| LEO 100 requests | 49.685 µs | 66.903 µs |
| GEO 24208 initialization | 776.77 ns | 749.64 ns |
| GEO single propagation | 670.05 ns | 620.77 ns |
| GEO 100 requests | 69.055 µs | 69.808 µs |
| Half-day 08195 initialization | 737.35 ns | 1072.4 ns |
| Half-day single propagation | 895.11 ns | 920.60 ns |
| Half-day 100 requests | 107.69 µs | 97.224 µs |

Some upstream measurements overlapped diagnostic compilation. LEO repeated
requests and half-day initialization have wide intervals (57.902–76.736 µs
and 956.31–1248.3 ns respectively). Treat these as preliminary measurements,
not reliable performance regression claims. Numerical correctness blocks
adoption independently of performance.

## License and validation status

The selected crate's exact license is MIT, copyright 2020 International Centre
for Neuromorphic Systems; its complete notice accompanies the diagnostic in
`LICENSE-sgp4-MIT`. It requires preserving the notice in redistributions.
Minimal backend dependencies are `chrono` (MIT OR Apache-2.0), `num-traits`
(MIT OR Apache-2.0), `libm` (MIT), and build dependency `autocfg`
(MIT OR Apache-2.0). Source/metadata were inspected in the registry packages.
No additional production dependency remains after the rejected experiment.
Siderust remains AGPL-3.0-only, version 0.15.0; any eventual commercial
redistribution must retain independent third-party notices. No upstream
implementation is attributed to Siderust or copied into its production source.

The native baseline and restored native code pass all nine numerical tests,
both pipeline tests, and both FFI tests (all three model codes and output-buffer
preservation on failure). The standalone audit passes Clippy. Workspace
formatting, diff whitespace checks, and
`cargo clippy --workspace --all-targets --all-features -- -D warnings` pass.
The full migration feature matrix, embedded checks, workspace tests/doctests,
and LLVM coverage were **not rerun as migration acceptance** after the
scientific failure. Prior PR coverage figures describe the native baseline,
not this rejected adapter. There is no new migration coverage result and no
claim that its definition of done has been met. All quality gates and fixtures
remain unchanged.
