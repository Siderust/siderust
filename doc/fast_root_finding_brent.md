# Fast Root Finding in Siderust (Brent + Endpoint Reuse + Native CPU)



## Background: what the benchmark is doing

The Criterion benchmark `find_night_periods_365day` ultimately calls:

- `calculus::solar::altitude_periods::find_sun_altitude_periods_via_culminations`

The heavy work happens in two phases:

1. **Culmination discovery** (upper/lower meridian crossings) using:
   - `calculus::events::find_dynamic_extremas::find_dynamic_extremas`

2. **Altitude boundary crossings** (twilight crossings) by bracketing crossings in each
   quasi-monotonic window and refining each root.

The dominant cost per function evaluation is computing the Sun’s apparent position,
which triggers a full **Earth VSOP87A evaluation** (thousands of trig terms), plus
nutation / precession corrections.

---

## The core functionality: derivative-free root refinement (Brent)

### Why bracketing methods matter here

Both culminating times and altitude crossings are obtained by **scanning** an interval
and detecting sign changes. That scan inherently provides **brackets**:

- You know the function is below the target at one endpoint and above at the other.

Once you have a bracket, you should prefer a **bracketing solver**.

### Newton vs Brent in this codebase

- Newton with finite differences needs ~3 expensive evaluations per iteration:
  - One evaluation at the current point
  - Plus two more at ±Δt to estimate the derivative

- Brent is derivative-free and typically needs ~1 new evaluation per iteration.

This is crucial when each evaluation is dominated by VSOP87 + trig.

### What was implemented

The root-finding module now contains separate implementations:

- `calculus::root_finding::newton`
- `calculus::root_finding::bisection`
- `calculus::root_finding::brent`

and exposes high-level APIs in `calculus::root_finding::mod`.

The main addition for performance-sensitive callers is:

- `brent::refine_root_with_values(...)`
- `root_finding::find_crossing_brent_with_values(...)`

These variants accept **pre-computed bracket endpoint values** so callers don’t need to
re-evaluate the expensive scalar function at the endpoints.

---

## Key optimization: reuse scan endpoint evaluations

### The problem

A typical event loop looks like:

1. Scan from `jd0` to `jd1`.
2. Evaluate `f(jd0)` and `f(jd1)` to detect a sign change.
3. If a sign change exists, call a root refiner.

A naive refiner will do this again:

- Evaluate `f(jd0)` and `f(jd1)` internally.

For expensive functions (Sun position), that’s **two redundant VSOP87 calls per event**.
Over a year horizon, those add up.

### The fix

Brent was extended with a variant that takes `f_lo` and `f_hi`:

- `refine_root_with_values(lo, hi, f_lo, f_hi, scalar_fn, threshold)`

The scan code already has `f_lo` and `f_hi` in hand (because it must compute them to
know whether a root is bracketed), so refinement starts without paying for those
evaluations again.

### Where it’s used

- Culmination refinement in `calculus::events::find_dynamic_extremas`
- Altitude crossing refinement in:
  - `calculus::events::altitude_periods`
  - `calculus::solar::altitude_periods`

---

## Culminations: scan with sin(), refine with hour-angle

Culmination discovery uses a coarse scan to locate brackets.

### Why `sin(H)` is used for detection

The hour angle $H$ is wrapped to $( -\pi, \pi ]$, so it has discontinuities.
Using `sin(H)` for sign-change detection is robust because it crosses zero at the
same events (upper/lower culminations), while avoiding wrap-related issues.

### Why refinement uses hour-angle directly

For refinement, solving on `H` (or `wrap_signed(H-π)` for lower culmination) is faster
than solving on `sin(H)`:

- Near the root, `sin(H)` flattens the slope slightly compared to `H`.
- Using `H` directly typically reduces iterations for the same time tolerance.

Implementation strategy:

- Detect brackets with `sin(H)` and `sin(H-π)`.
- Refine roots with Brent on:
  - `H(jd)` for upper culmination
  - `wrap_signed(H(jd) - π)` for lower culmination

while passing precomputed endpoint values into `refine_root_with_values`.

---

## Accuracy vs speed: tolerance selection

Brent’s termination uses a time-domain tolerance on the bracket size.

A very small tolerance (e.g., `1e-11` days) can be **overkill** for this workload and
forces extra iterations, increasing VSOP87 calls.

The tuned tolerance is:

- `TOLERANCE = 1e-9 days` (≈ 86 μs)

This is:

- Tight enough to satisfy the existing unit/integration tests for root location.
- Loose enough to avoid wasting iterations in event-finding loops.

If you need stricter event timing (e.g., sub-10 μs), expect a measurable runtime
increase.

---

## Why it is so fast (the real reason)

The performance wins come from two independent effects:

### 1) Fewer expensive evaluations

- Replacing Newton-with-finite-diff with Brent removes the ±Δt derivative probes.
- Reusing scan endpoint values avoids **two wasted evaluations per refined root**.

Because the scalar functions are dominated by VSOP87A trig, reducing the number of
calls is the highest-leverage optimization.

### 2) Better CPU codegen with native target features

The hot loops are dominated by vectorized trig and multiply-add accumulation.

On modern CPUs (e.g., Ryzen with AVX2/FMA), compiling with native CPU features can:

- Enable true FMA lowering for `mul_add`.
- Improve SIMD code paths.
- Reduce libcall/helper overhead that shows up in profilers.

For local benchmarking and deployments where portability is not required, use:

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench solar_altitude find_night_periods_365day
```

In measurements on this workspace:

- Default build: ~1.15 s
- With `target-cpu=native`: ~0.55 s

(Exact numbers depend on CPU, Rust version, and background load.)

---

## How to reproduce

### Baseline benchmark

```bash
cargo bench --bench solar_altitude find_night_periods_365day
```

### Fast benchmark (native CPU)

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench solar_altitude find_night_periods_365day
```

### Notes

- Criterion defaults aim for stable statistics, not minimal wall time. If you only
  want a quick datapoint, you can reduce the sample count via Criterion config or
  just run once under `cargo run --release` with a dedicated timing harness.

---

## Trade-offs and safety notes

- Bracketing methods (Brent/Bisection) require valid brackets. The scan code must
  be conservative enough to not miss sign changes.
- The endpoint-value reuse assumes that `f_lo` and `f_hi` correspond to
  `scalar_fn(lo) - threshold` and `scalar_fn(hi) - threshold` exactly.
- `target-cpu=native` produces binaries optimized for the build machine.
  Do not ship these binaries to unknown CPUs.

---

## Related files

- Root finding:
  - `src/calculus/root_finding/brent.rs`
  - `src/calculus/root_finding/mod.rs`
- Culminations:
  - `src/calculus/events/find_dynamic_extremas.rs`
- Altitude crossings:
  - `src/calculus/events/altitude_periods.rs`
  - `src/calculus/solar/altitude_periods.rs`
- Perf context:
  - `doc/solar_altitude_perf_report.md`
