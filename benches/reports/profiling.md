# Profile a benchmark on Linux

First reproduce the timing difference with Criterion on a controlled host.
Select one case and keep its input, features and dataset fixed. Criterion's
`--profile-time` runs the selected workload repeatedly without collecting
the normal statistical estimates.

## CPU stacks and counters

Build with line tables and use the executable path printed by Cargo.

```sh
CARGO_PROFILE_BENCH_DEBUG=line-tables-only \
  cargo bench --locked -p siderust --bench solar_altitude --features bench-internals --no-run
```

Set `bench_bin` to the executable path printed above. These commands scope
profiling to the process and use userspace events.

```sh
case_id='public_api/solar/30d/below_threshold/astro_night'
perf stat -e task-clock,cycles:u,instructions:u -- \
  "$bench_bin" --bench --profile-time 10 --exact "$case_id"
perf record -e cycles:u --call-graph dwarf -- \
  "$bench_bin" --bench --profile-time 10 --exact "$case_id"
perf report
```

Use counters to explore instruction cost and CPU activity. Use sampled stacks
to find where the workload spends CPU time. Counters cover the entire profiled
process, including startup and harness work. DWARF unwinding can increase
sampling overhead. If counters are restricted or unavailable, try a permitted
software event such as `cpu-clock:u` for stacks or use an approved profiling
host. See the [perf record manual](https://man7.org/linux/man-pages/man1/perf-record.1.html).

With [cargo-flamegraph](https://github.com/flamegraph-rs/flamegraph) installed,
the equivalent filtered workflow is:

```sh
CARGO_PROFILE_BENCH_DEBUG=line-tables-only \
  cargo flamegraph -p siderust --bench solar_altitude --features bench-internals -- \
  --bench --profile-time 10 --exact 'public_api/solar/30d/below_threshold/astro_night'
```

Frame width represents sample weight including callees. Inspect both the
top stacks and the callers leading to them. Re-run Criterion without
instrumentation to measure the impact of a change.

## Heap use

Launch [Heaptrack](https://github.com/KDE/heaptrack) on the same symbolized
executable. Use the trace path it reports to inspect allocation sites, volume
and peak live heap.

```sh
heaptrack "$bench_bin" --bench --profile-time 10 --exact "$case_id"
```

Compare identical configurations and account for Criterion buffers and
process initialization. Heap profiles do not quantify flash or stack use.
For resident memory use an OS process metric with its accounting definition.
For device memory use a section map and measure stack and heap high-water
on the board. The [compact footprint report](compact_ephemeris.md) provides
an example of build-size accounting.
