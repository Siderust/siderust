# Solar Altitude bench — perf investigation (`find_night_periods_365day`)

## Scope / inputs

- **Benchmark file**: `benches/solar_altitude.rs`
- **Profiled benchmark**: `find_night_periods_365day` (Criterion filter argument)
- **Profiler command** (from `perf.data` header):
  - `perf record -F 997 --call-graph dwarf,64000 -g -o perf.data …/solar_altitude-… find_night_periods_365day`
- **Capture date**: Thu **2026-02-05**
- **CPU**: AMD Ryzen 7 5700X (AVX2/FMA capable; perf shows SSE2 path in hot code)

Goal: identify where the benchmark spends CPU cycles and list the highest-leverage ways to reduce runtime (no code changes in this pass).

---

## Executive summary (what’s “wasting time”)

Perf shows the benchmark is dominated by **ephemeris + trig**, specifically:

1. **~74% of cycles** in `siderust::calculus::vsop87::vsop87_impl::position` (Earth VSOP87A evaluation via coordinate transforms).
2. Inside that, most time is `wide::f64x4::cos/sin_cos` (range reduction + polynomial trig).
3. A secondary but visible cost is **nutation** trig (`libm::__sincos_fma` + `astro::nutation::get_nutation`).
4. **~9% of cycles** in `fma` helper symbols suggests `mul_add` is not compiling down to a single FMA instruction in this build (or is going through a dispatch/libcall path).

Put simply: the bench does *a lot* of Sun position evaluations, and each Sun position requires an expensive Earth VSOP87A series evaluation with thousands of cosines.

---

## Top costs from perf (cycles)

From `perf report` (2K samples):

| Share (approx.) | Symbol | Meaning |
|---:|---|---|
| ~74% | `siderust::calculus::vsop87::vsop87_impl::position` | VSOP87 position loop (called from `Earth::vsop87a`) |
| ~64% (within VSOP) | `wide::f64x4::cos` / `wide::f64x4::sin_cos` (inlined) | SIMD trig dominates the VSOP loop |
| ~10% | `libm.so.6::__sincos_fma` | scalar trig (mostly nutation) |
| ~6.5% + ~2.6% | `fma` + `compiler_builtins … fma_with_fma` | `mul_add` fused semantics implemented via helper |
| ~2–7% | `siderust::astro::nutation::get_nutation` | nutation series / argument setup (plus its trig) |

This is a compute-bound profile; the hotspots are math kernels, not allocations or IO.

---

## How `find_night_periods_365day` reaches VSOP87 (call-path explanation)

### Benchmark entrypoint

`benches/solar_altitude.rs` benchmarks `find_night_periods(site, period, twilight)` which calls:

- `src/calculus/solar/altitude_periods.rs::find_sun_altitude_periods_via_culminations`

### Two heavy phases inside `find_sun_altitude_periods_via_culminations`

1. **Culmination discovery** (`find_dynamic_extremas`):
   - Calls a closure `get_equatorial(jd)` repeatedly.
   - That closure starts from a heliocentric “Sun at origin” position and does:
     - `helio = cartesian::position::Ecliptic::<AU>::CENTER`
     - `helio.transform(jd)` → geocentric equatorial mean J2000 cartesian
     - cartesian→spherical, then precession + nutation RA correction

2. **Altitude threshold crossings** (sun below twilight):
   - Repeatedly evaluates `altitude_fn(jd) = sun_altitude_rad(jd, &site)`
   - `sun_altitude_rad` calls `Sun::get_horizontal(jd, site)` (`src/calculus/solar/sun_equations.rs`),
     which again uses `helio.transform(jd)` as its first step.

### Critical detail: `helio.transform(jd)` triggers `Earth::vsop87a(jd)`

The translation from **Heliocentric → Geocentric** for positions is implemented in:

- `src/coordinates/transform/centers/position/to_geocentric.rs`

For `Position<Heliocentric, …>` it does:

- `earth_helio_ecl_au = Earth::vsop87a(jd).position`
- rotate Earth into the requested frame
- translate: `geocentric = heliocentric - earth_position`

So every time code transforms the heliocentric origin to geocentric (i.e., “compute Sun’s geocentric direction”), it pays for one **full Earth VSOP87A** evaluation.

---

## Why VSOP87 dominates (and what inside it dominates)

### 1) Earth VSOP87A is thousands of cosine terms per call

From the generated VSOP87A dataset in this workspace (term counts visible in `target/release/build/*/out/vsop87a.rs`):

- `EARTH_X*` total terms: **1,577**
- `EARTH_Y*` total terms: **1,590**
- `EARTH_Z*` total terms: **371**
- **Total**: **3,538** VSOP87 terms per Earth position evaluation

Each VSOP87 term is `a * cos(b + c*T)`, so that’s ~3,538 cosine evaluations per `Earth::vsop87a(jd)`.

Repro command for the term counts:

```bash
rg -n "pub static EARTH_[XYZ][0-5]: \\[Vsop87;" target/release/build/*/out/vsop87a.rs
```

### 2) The algorithm calls Sun position *a lot* for a 365-day horizon

At a high level, the 365-day case implies:

- A coarse scan for culminations (20-minute step in `find_dynamic_extremas`) → ~72 steps/day → ~26k evaluations.
- Newton refinements for each culmination and each altitude crossing.
  - Both Newton implementations use central finite differences (±1s), which multiplies expensive evaluations by ~3 per iteration.

That makes the total number of VSOP87 evaluations easily reach **tens of thousands** for a year.

Even conservative back-of-envelope:

> 50k ephemeris evaluations × 3,538 terms ≈ **177 million cosines**

This aligns with perf: the VSOP87 position kernel dominates.

### 3) SIMD trig cost is the bottleneck inside VSOP87

Within `vsop87_impl::coord_value`, perf shows most time inside `wide::f64x4::cos/sin_cos`.
The inlined leaf frames (`blend`, `round_int`, shifts, xor/and) are characteristic of:

- argument range reduction, and
- polynomial approximation.

This is expected and essentially means “trig is expensive”.

### 4) `mul_add` lowering shows up as `fma()` helper cost

You intentionally use `mul_add` in the VSOP87 loop (argument formation and series accumulation).
Perf showing `fma` helper symbols at ~9% suggests this build is not fully taking advantage of a compile-time `+fma` target feature.

This is an actionable lever for benchmarks and non-portable builds.

---

## Improvement options (ordered by impact)

### 1) Reduce the *number* of expensive evaluations (highest leverage)

**A. Replace Newton+finite-difference with a derivative-free bracketing solver (Brent)**

Both:

- `src/calculus/root_finding.rs::refine_root_newton` (threshold crossings), and
- `src/calculus/events/find_dynamic_extremas.rs` refinement loop (culminations)

spend **3 expensive evaluations per iteration** due to ±1s finite-difference derivatives.

But you already have brackets from scanning sign changes, so you can use **Brent’s method**:

- robust like bisection,
- faster convergence in practice,
- typically ~1 new function evaluation per iteration,
- no ±1s side evaluations.

Impact: large, because it directly reduces VSOP87 calls.

**B. Sun-specific shortcut for culminations**

`find_dynamic_extremas` is generic and scans 365 days at 20-minute granularity.
For the Sun, upper culmination (solar noon) can be approximated from:

- longitude, and
- equation of time (EoT),

then locally refined.

Impact: can remove most of the coarse scanning cost for the Sun case.

**C. Two-tier model: fast scan + accurate refine**

For twilight windows you usually need accurate event *times*, not high-accuracy positions at every intermediate step.

A standard strategy:

1. Use a **fast solar model** (few trig terms) to scan/bracket crossings.
2. Refine final roots using the full model (VSOP87 + nutation + parallax) only near candidate roots.

Impact: potentially orders of magnitude, depending on how aggressive the “fast” model is.

### 2) Make each VSOP87 evaluation cheaper (still valuable)

**A. Compile benches with native CPU features (`target-cpu=native`)**

Perf shows SSE2 intrinsics (`_mm_*`) and non-trivial time in `fma` helper symbols.
On Ryzen 5700X, AVX2 and FMA are available.

For benchmarking / local deployment, try:

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench solar_altitude -- find_night_periods_365day
```

Expected effects:

- reduced/eliminated `fma` helper overhead,
- potentially wider SIMD and faster `wide` operations,
- lower cycles/term in VSOP87.

**B. Reduce VSOP87 term count for this use-case**

If your application tolerates small angular error in the Sun’s position (common for sunrise/sunset):

- truncate the series (drop higher powers and/or smallest-amplitude terms), or
- switch to a lower-order solar/earth model.

This is the main path to “big-O” savings *per evaluation*.

**C. Reduce `coord_value` SIMD overhead (few‑percent class)**

Perf shows some overhead in:

- `bytemuck::cast` / `read_unaligned` inside `wide::to_array()`.

Avoiding `to_array()` and doing a vector dot + horizontal sum would reclaim some cycles.
This won’t beat “reduce trig calls”, but it’s low-hanging once the algorithmic side is addressed.

### 3) Secondary: reduce nutation work

Nutation isn’t the primary bottleneck, but it is visible (~10% scalar trig + nutation frames).

Options:

- cache nutation results (or rotation matrices) at a coarser cadence (hourly/daily) and interpolate,
- use a shorter nutation model during scan steps.

---

## Suggested next measurements (to decide what to implement)

1. **Measure raw VSOP cost**:
   - `cargo bench --bench vsop87` (already exists) to quantify `Earth::vsop87a` baseline.
2. **Quantify “native CPU features” impact**:
   - compare default vs `RUSTFLAGS="-C target-cpu=native"` for the same benchmark.
3. **Estimate evaluation counts**:
   - count how many times `sun_altitude_rad` is evaluated for the 365‑day case (instrumentation or log counters in a profiling build).
4. **Re-profile after any change**:
   - success criterion: shrink `vsop87_impl::position` width first; everything else is secondary.
