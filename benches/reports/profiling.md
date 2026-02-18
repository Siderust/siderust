# Profiling Benchmarks (Criterion + `perf` + Flamegraphs)

> Note: This report lives under `benches/reports/` because it documents
> benchmark-driven performance work.

This document is a practical workflow for turning a benchmark regression into
actionable profiling data. The intended loop is:

1) measure with Criterion,  
2) explain with a flamegraph (or `perf report`),  
3) optimize,  
4) re-measure.

## What flamegraphs do (and do not) tell you

A flamegraph is a visualization of sampled CPU stacks. The only dimension that
correlates with “time spent” is **width**: a wider frame means the profiler saw
that function (inclusive of callees) more often.

Flamegraphs are great for answering “which code dominates CPU time?” They are
not a wall-clock timing tool, and they do not automatically explain allocation
or IO behavior. Treat them as a map that tells you where to look next.

## Prerequisites (Linux)

Install the Rust wrapper:

```bash
cargo install flamegraph
```

You also need `perf` available on the system. On Debian/Ubuntu:

```bash
sudo apt install linux-tools-common linux-tools-$(uname -r)
```

On Arch:

```bash
sudo pacman -S perf
```

## Kernel permissions for `perf`

Many distros restrict performance counters by default. If `perf` fails with
errors mentioning `perf_event_paranoid`, relax the settings for development.

Temporary (until reboot):

```bash
sudo sysctl -w kernel.perf_event_paranoid=1
sudo sysctl -w kernel.kptr_restrict=0
```

Permanent (optional): create `/etc/sysctl.d/99-perf.conf` with:

```conf
kernel.perf_event_paranoid = 1
kernel.kptr_restrict = 0
```

and apply:

```bash
sudo sysctl --system
```

If you cannot change sysctl globally, another option is to grant capabilities:

```bash
sudo setcap cap_perfmon,cap_sys_ptrace=eip $(which perf)
```

## Build configuration: keep symbols in release

Profiling optimized code without symbols is usually a waste of time. Keep DWARF
symbols in release builds so stacks resolve to function names.

Add to `Cargo.toml`:

```toml
[profile.release]
debug = true
```

This keeps optimizations enabled; it just includes enough metadata to make
profilers useful.

## Profiling a Criterion benchmark

To profile a whole benchmark binary:

```bash
cargo flamegraph --bench solar_altitude
```

Criterion supports filtering. Arguments after `--` are passed through to the
benchmark harness:

```bash
cargo flamegraph --bench moon_altitude -- --bench find_moon_above_horizon_365day
```

If you only care about one case (e.g. the 365‑day horizon), filtering is the
single best way to reduce noise and iteration time.

## Reducing noise from Criterion itself

Criterion’s default configuration aims for stable statistics, which means it
may run more iterations than you want while profiling.

For deep profiling sessions, it can be worth temporarily lowering `sample_size`
and warm-up duration in the benchmark file, or creating a dedicated micro-bench
that calls exactly the code path you want to study.

When you do that, make sure you still prevent dead-code elimination:

```rust
use criterion::black_box;
black_box(your_function(black_box(input)));
```

## Common pitfalls

- Profiling `debug` builds (results rarely generalize to release performance).
- Reading flamegraph “height” as “time” (width is what matters).
- Optimizing thin frames before removing the dominant hotspot.
- Making a change without re-running Criterion to confirm impact.

