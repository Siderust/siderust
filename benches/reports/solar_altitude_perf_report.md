# Solar Altitude Bench — Architecture and Performance Notes

> Documents benchmark-driven performance work for `benches/solar_altitude.rs`.

## What is measured

The benchmark groups compare solar threshold event discovery over 30-, 184-, and 365-day
windows at Roque de los Muchachos:

- `solar/daily_predictor/*` — public Option A API (`below_threshold`, `altitude_ranges`, `crossings`)
- `solar/chebyshev_fallback_baseline/*` — internal Chebyshev-first engine (daily predictor disabled)
- `solar/scan_brent_baseline/*` — internal scan+Brent baseline via `bench_internals`

Thresholds exercised in the daily-predictor group include astronomical twilight (−18°),
nautical bands, and horizon crossings.

## Engine architecture (current branch)

### Default: solar daily predictor

1. Iterate full MJD days overlapping the requested window.
2. Build an **analytic daily state** (approximate RA/Dec/transit) from a compact solar model.
3. Predict morning/evening crossing candidates from hour-angle geometry.
4. Expand local brackets (15 min → 4 h) and refine with **`sun_altitude_rad` + Brent**.
5. Accept a day only when every in-window candidate refines successfully.

Performance depends primarily on **precise `sun_altitude_rad` evaluation count**,
not on scanning the full window at a fixed step.

### Fallback: local Chebyshev-first

When the daily model is unreliable (grazing latitude, partial-window edge cases, bracket
failures, or failed always-above/below safety checks), the clipped day falls back to the
generic Chebyshev-first crossing engine in `src/event/search/crossings.rs`.

### Final fallback / baseline: local scan+Brent

Per-segment or whole-window scan+Brent is used only as an internal baseline (tests/benches)
or when Chebyshev marks a segment unsafe.

## Expected diagnostics (normal mid-latitude cases)

For ordinary sites and twilight thresholds:

- `scan_fallback_days == 0`
- `bracket_failures == 0` (or handled by local day fallback without silent partial results)
- `precise_evaluations` substantially lower than Chebyshev-first or scan+Brent baselines

Polar and grazing cases should show non-zero `chebyshev_fallback_days` and still match
scan baselines in tests.

## How to run

```bash
cargo bench --features bench-internals --bench solar_altitude
cargo bench --no-run --features bench-internals
```

Plain `cargo bench --no-run` skips these benches because they require `bench-internals`.

## Speedup claims

This report does **not** claim measured speedups unless benchmark numbers are pasted here
after a local run. Compare groups `solar/daily_predictor` vs `solar/chebyshev_fallback_baseline`
vs `solar/scan_brent_baseline` for the same window/threshold.
