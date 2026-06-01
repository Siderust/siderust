# Benchmarks

This crate uses [Criterion](https://crates.io/crates/criterion) for benchmarks (stable Rust).
Results go to `target/criterion/` with HTML reports.

## Quick Start

```bash
cargo bench                           # run Criterion suites (no BSP required)
SIDERUST_BSP_PATH=/path/to/de440.bsp cargo bench --bench de441   # runtime JPL
```

## Comparative Suites

These run the **same operations** across different backends or body types so
Criterion's HTML reports show relative performance in a single chart.

| Benchmark | What it compares | Command |
|-----------|-----------------|---------|
| `ephemeris_comparison` | VSOP87 `Ephemeris` trait methods | `cargo bench --bench ephemeris_comparison` |
| `altitude_comparison`  | Sun vs Moon vs Star (Sirius) for single-point, 7-day, 30-day, 365-day searches | `cargo bench --bench altitude_comparison` |

## Per-Module Benchmarks

Detailed benchmarks for each subsystem with multiple time horizons and algorithm variants.

| Benchmark | Scope | Command |
|-----------|-------|---------|
| `converters` | Cartesian ↔ Spherical conversion | `cargo bench --bench converters` |
| `vsop87` | VSOP87A/E evaluation for all 8 planets | `cargo bench --bench vsop87` |
| `elp2000` | ELP2000-82B Moon position at multiple epochs | `cargo bench --bench elp2000` |
| `solar_altitude` | Sun night-period finding (1/7/30/365-day) | `cargo bench --bench solar_altitude` |
| `moon_altitude` | Moon above/below horizon, altitude ranges, algorithm comparison | `cargo bench --bench moon_altitude` |
| `star_altitude` | Fixed-star altitude, thresholds, crossings | `cargo bench --bench star_altitude` |
| `de441` | Runtime DE4xx ephemeris body queries (needs `SIDERUST_BSP_PATH`) | `cargo bench --bench de441` |

## Filtering

Criterion accepts a filter string (substring/regex) after the bench name:

```bash
cargo bench --bench moon_altitude compute_altitude
cargo bench --bench ephemeris_comparison -- --help
cargo bench --bench ephemeris_comparison -- --sample-size 50
```

## Reports

Benchmark-driven performance investigations and profiling notes live under:

- `benches/reports/README.md`
