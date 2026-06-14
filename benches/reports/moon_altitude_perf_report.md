# Moon Altitude Bench — Architecture and Performance Notes

> Documents benchmark-driven performance work for `benches/moon_altitude.rs`.

## Benchmark groups

```
public_api/moon/{30d,184d,365d}
  above_threshold/horizon              — above_threshold(&Moon, …, 0°)
  below_threshold/horizon              — below_threshold(&Moon, …, 0°)
  altitude_ranges/observation_window   — altitude_ranges(&Moon, …, 0°, 30°)

engines/moon/{30d,184d,365d}
  above_threshold/chebyshev_context    — default engine (same as public_api)
  above_threshold/scan_brent_baseline  — uniform scan + Brent (internal baseline)
  altitude_ranges/chebyshev_context    — default engine (same as public_api)
  altitude_ranges/scan_brent_baseline  — uniform scan + Brent (internal baseline)
```

All groups require `bench-internals`. Run with:

```bash
cargo bench --features bench-internals --bench moon_altitude
cargo bench --no-run --features bench-internals
```

Plain `cargo bench --no-run` skips these benches because they require `bench-internals`.

## Engine architecture

### Default: `MoonAltitudeContext` + Chebyshev-first crossings

1. Build a cached lunar context for the search window (ELP2000 + nutation caches).
2. Fit Chebyshev segments on `sin(altitude) − sin(threshold)`.
3. Validate roots against the precise lunar altitude signal.
4. Fall back per unsafe segment to local scan+Brent.

Topocentric parallax and Earth rotation stay in the precise evaluator. The Chebyshev layer
reduces the number of expensive evaluations needed across long windows.

### Scan+Brent baseline (internal)

Uniform scan + Brent refinement over the full window. Not exposed in the public API.

## Public API

All user-facing calls use `SearchOpts::default()`:

```rust
above_threshold(&Moon, &site, window, threshold, SearchOpts::default())
below_threshold(&Moon, &site, window, threshold, SearchOpts::default())
altitude_ranges(&Moon, &site, window, h_min, h_max, SearchOpts::default())
```

No scan-step or algorithm knobs are available publicly.

## Speedup claims

This report does **not** claim measured speedups unless benchmark numbers are pasted here
after a local run. Compare `engines/moon/{N}/above_threshold/chebyshev_context` vs
`scan_brent_baseline` for the same window to measure the speedup.
