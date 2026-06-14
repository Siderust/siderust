# Solar Altitude Bench — Architecture and Performance Notes

> Documents benchmark-driven performance work for `benches/solar_altitude.rs`.

## Benchmark groups

```
public_api/solar/{30d,184d,365d}
  below_threshold/astro_night          — below_threshold(&Sun, …, −18°)
  above_threshold/horizon              — above_threshold(&Sun, …, 0°)
  altitude_ranges/astro_to_nautical    — altitude_ranges(&Sun, …, −18°, −12°)  [batch daily path]
  crossings/horizon                    — crossings(&Sun, …, 0°)

engines/solar/{30d,184d,365d}
  below_threshold/daily_predictor      — default engine (same as public_api)
  below_threshold/chebyshev_baseline   — Chebyshev-first, daily predictor disabled
  below_threshold/scan_brent_baseline  — uniform scan + Brent (internal baseline)
  altitude_ranges/daily_batch          — batch daily path (same as public_api)
  altitude_ranges/chebyshev_baseline   — Chebyshev-first, no daily predictor
  altitude_ranges/scan_brent_baseline  — uniform scan + Brent (internal baseline)
```

All groups require `bench-internals`. Run with:

```bash
cargo bench --features bench-internals --bench solar_altitude
cargo bench --no-run --features bench-internals
```

Plain `cargo bench --no-run` skips these benches because they require `bench-internals`.

## Engine architecture

### Default: solar daily predictor

1. Iterate full MJD days overlapping the requested window.
2. Build an **analytic daily state** (approximate RA/Dec/transit) from a compact solar model.
3. Predict morning/evening crossing candidates from hour-angle geometry.
4. **Fast-model bracket** (Meeus-style `sin(altitude)` — ~8 arcmin accuracy, ~10× cheaper than
   the precise model) to narrow the root location.
5. **Newton/secant polish** (up to 6 iterations) using precise `sun_altitude_rad` residuals,
   starting from the fast root. Accepts when `|residual| ≤ 1e-9` (time accuracy ~16–30 µs,
   well within the declared 86 µs `time_tolerance`).
6. Fallback to full-model Brent if Newton/secant does not converge (polar/grazing cases).

For `altitude_ranges`, a **batch daily path** computes both threshold crossings (h_min and h_max)
in a single day loop, sharing the daily state computation.

Performance scales with **precise `sun_altitude_rad` evaluation count**, not window size.

### Chebyshev-first baseline (internal)

Generic Chebyshev polynomial fit on `sin(altitude) − sin(threshold)` per segment,
with per-segment scan+Brent fallback. Daily predictor disabled.

### Scan+Brent baseline (internal)

Uniform scan at 2-hour steps across the full window, Brent refinement at each sign change.
Not exposed in the public API.

## Expected diagnostics (normal mid-latitude cases)

For ordinary sites and twilight thresholds (measured: Roque de los Muchachos, 30 days, −18°):

| Metric | Value |
|--------|-------|
| `crossings` | 60 |
| `precise_evaluations` | ~261 |
| `newton_accepted` | ~51 / 60 |
| `evals/crossing` | ~4.3 |
| `bracket_failures` | 0 |
| `scan_fallback_days` | 0 |

Polar and grazing cases will show non-zero `chebyshev_fallback_days` and still match
scan baselines in tests.

## Measured performance (Roque de los Muchachos, 2026, release build)

Measured on Linux 6.17, AMD Ryzen, `cargo bench --features bench-internals`.

| Benchmark | Time | Notes |
|-----------|------|-------|
| `public_api/solar/365d/below_threshold/astro_night` | **471 ms** | −45% vs pre-optimization 0.85 s |
| `public_api/solar/365d/above_threshold/horizon` | **467 ms** | |
| `public_api/solar/365d/altitude_ranges/astro_to_nautical` | **955 ms** | −40% vs pre-optimization 1.58 s |
| `public_api/solar/365d/crossings/horizon` | **552 ms** | |

The speedup comes from reducing precise `sun_altitude_rad` evaluations per root from ~9 (old
Brent-only path) to ~4.3 (Newton/secant polish accepts 85% of mid-latitude crossings in 3
iterations). No public API changes; no reduction in crossing accuracy.
