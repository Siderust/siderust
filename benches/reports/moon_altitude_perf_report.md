# Moon Altitude Bench — Architecture and Performance Notes

> Documents benchmark-driven performance work for `benches/moon_altitude.rs`.

## What is measured

The benchmark groups compare lunar above/below-horizon window finding over 30-, 184-, and
365-day windows at Roque de los Muchachos:

- `moon/chebyshev_context/*` — public Option A API (`above_threshold`, `below_threshold`)
- `moon/scan_brent_baseline/*` — internal scan+Brent baseline via `bench_internals`

## Engine architecture (current branch)

### Default: `MoonAltitudeContext` + Chebyshev-first crossings

1. Build a cached lunar context for the search window (ELP2000 + nutation caches).
2. Fit Chebyshev segments on `sin(altitude) − sin(threshold)`.
3. Validate roots against the precise lunar altitude signal.
4. Fall back per unsafe segment to local scan+Brent.

Topocentric parallax and Earth rotation remain in the precise evaluator; the Chebyshev layer
reduces the number of expensive samples needed across long windows.

### Fallback / baseline: scan+Brent

Internal baseline forces uniform scan+Brent over the window for comparison in tests and benches.
This is not exposed in the public API (`SearchOpts` contains only `time_tolerance`).

## Public API

All user-facing calls use:

```rust
above_threshold(&Moon, &site, window, threshold, SearchOpts::default())
below_threshold(&Moon, &site, window, threshold, SearchOpts::default())
```

No scan-step or algorithm knobs are available publicly.

## How to run

```bash
cargo bench --features bench-internals --bench moon_altitude
cargo bench --no-run --features bench-internals
```

## Speedup claims

This report does **not** claim measured speedups unless benchmark numbers are pasted here
after a local run.
