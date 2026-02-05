# Benchmarks

This crate uses [Criterion](https://crates.io/crates/criterion) for benchmarks (so it works on stable Rust).

## Run all benches

From the repository root:

```bash
cargo bench
```

Results are written under `target/criterion/`.

## Run a single bench binary

```bash
cargo bench --bench vsop87
cargo bench --bench converters
cargo bench --bench solar_altitude
cargo bench --bench moon_altitude
cargo bench --bench elp2000
```

## Filter to specific benchmarks inside a bench

Criterion accepts an optional filter string (substring/regex) after the bench name:

```bash
cargo bench --bench moon_altitude compute_altitude
```

To pass Criterion CLI options, add `--`:

```bash
cargo bench --bench moon_altitude -- --help
cargo bench --bench moon_altitude -- --sample-size 20
```

## What’s available

- `converters`: cartesian/spherical conversion benchmarks.
- `vsop87`: VSOP87 series evaluation benchmarks.
- `solar_altitude`: solar altitude period-finding benchmarks.
- `moon_altitude`: lunar altitude + horizon/altitude search benchmarks.
- `elp2000`: ELP2000 lunar position evaluation benchmarks.

