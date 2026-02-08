# Moon altitude bench — perf investigation (`find_moon_above_horizon_365day`)

## Scope / inputs

- **Benchmark file**: `benches/moon_altitude.rs`
- **Profiled benchmark**: `find_moon_above_horizon_365day` (Criterion filter argument)
- **Capture date**: Thu **2026-02-05**
- **Run command** (local):

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench moon_altitude
```

Goal: get Moon altitude *period finding* over **365 days** as fast as practical while keeping physically-correct topocentric altitude (parallax matters for the Moon).

---

## Executive summary

This benchmark is dominated by **ELP2000 Moon ephemeris evaluations** plus coordinate transforms needed to produce **topocentric horizontal altitude**.

Key improvements came from reducing the **number of expensive altitude evaluations** (not from micro-optimizing math):

- Replace the generic high-resolution scan strategy with a Moon-specific coarse scan.
- Refine only bracketed crossings.
- Avoid extra altitude calls during crossing classification.

Current results (Ryzen-class CPU, `target-cpu=native`):

- `moon_above_horizon/find_moon_above_horizon_365day`: **~1.398 s**
- `moon_algorithm_comparison/moon_above_horizon_scan_365day` (generic baseline): **~13.66 s**
- `moon_altitude_single/compute_altitude`: **~241 µs**

So the optimized path is ~**9.8× faster** than the generic scan baseline.

---

## What “expensive” means for the Moon

A single altitude evaluation uses:

1. `Moon::get_geo_position(jd)` (ELP2000-82B; the dominant cost)
2. Frame transforms: ecliptic → equatorial J2000
3. **Topocentric correction** (`to_topocentric(site, jd)`) for lunar parallax
4. Precession J2000 → mean-of-date
5. Equatorial → horizontal for the observer
6. Cartesian → spherical to extract `alt`

This is implemented in `src/calculus/lunar/altitude_periods.rs::moon_altitude_rad`.

The Moon must be handled differently than the Sun:

- Solar code relies on VSOP87 (fast-ish per call but many trig terms).
- Lunar code uses ELP2000 (heavier per call), and the Moon’s strong parallax makes topocentric calculations non-optional.

---

## Algorithm overview (final)

The optimized public wrappers:

- `find_moon_above_horizon(site, period, threshold)`
- `find_moon_below_horizon(site, period, threshold)`

use:

- `find_moon_altitude_periods_fast(site, period, condition)`

### 1) Coarse scan (2-hour step)

We scan with a fixed step (`MOON_ALTITUDE_SCAN_STEP = 2 hours`) and detect **sign changes** in:

- `f(jd) = altitude(jd) - threshold`

Each sign change brackets a crossing (rise/set or threshold crossing). This collapses work from “evaluate every few minutes” to “evaluate ~12 times/day”.

### 2) Root refinement only when bracketed

Each bracketed crossing is refined with a local Brent-style solver.

Important detail: refinement reuses the already-computed endpoint values `f(lo)` and `f(hi)`.

### 3) Crossing classification without extra altitude calls

The original direct-scan implementation classified each crossing by evaluating altitude before/after the root (two extra ephemeris calls per crossing).

The final approach labels each crossing using the **scan-time sign direction**:

- If `prev_f < 0` and `next_f > 0`, altitude is increasing through the boundary → “rising through threshold”.

For `Above(threshold)` this corresponds directly to **entering** the valid region; for `Below(threshold)` the direction is inverted.

This avoids ~2 additional `moon_altitude_rad` calls per crossing and was a major speed win.

### 4) Relaxed time tolerance for lunar rise/set

For long-range planning, sub-millisecond root tolerance is unnecessary. The lunar solver uses a relaxed tolerance of ~**2 minutes** to cut refinement iterations while still matching the scan/cross-check tests.

---

## Benchmark results

From the last run on 2026-02-05:

- **Single altitude** (`moon_altitude_single/compute_altitude`): ~**241 µs**
- **1 day** above horizon: ~**4.10 ms**
- **7 day** above horizon: ~**27.3 ms**
- **30 day** above horizon: ~**115.8 ms**
- **365 day** above horizon: ~**1.398 s**

Algorithm comparison:

- **Generic scan** (`moon_algorithm_comparison/moon_above_horizon_scan_365day`): ~**13.66 s**
- **Optimized** (`moon_algorithm_comparison/moon_above_horizon_culminations_365day`): ~**1.40 s**

Note: the “culminations” benchmark name reflects the historical implementation. Current optimized behavior is driven by the fast scan + refined crossings path.

---

## Why < 1 second is hard here

Even with a coarse scan, the year-long search still needs to find and refine hundreds of rise/set crossings.

Given an altitude evaluation cost of ~241 µs:

- The theoretical lower bound is on the order of **(number of altitude calls) × 241 µs**.

Getting under 1 second likely requires either:

1. A **two-tier ephemeris** (cheap approximate model for bracketing; ELP2000 only for final refinement), or
2. Parallelism (not currently used; would require API changes or internal threading), or
3. Significant simplification/caching inside the ELP2000 evaluation + transform pipeline.

---

## Suggested next steps (if we want < 1s)

Highest-leverage ideas:

1. **Two-tier bracketing**:
   - Use a cheap lunar approximation (few terms) to bracket crossings.
   - Refine final roots with full `moon_altitude_rad`.

2. **Cache slow-changing transforms**:
   - Precession/nutation matrices and/or sidereal time terms can be cached per hour/day and interpolated.

3. **Reduce work per refinement**:
   - Keep tolerance at ~2 minutes for planning-grade outputs.
   - Consider bisection with fixed iteration cap when brackets are tight.

4. **Instrument evaluation counts**:
   - Add counters (behind a feature flag) for number of altitude evaluations per year search to confirm where time goes.
