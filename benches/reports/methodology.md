# Measurement and reproducibility

## Define the operation

Specify the public call or internal path, work per iteration, output and units.
For an ephemeris query record body, origin, frame, epoch and time scale. For an
event search also record site, window, thresholds, search options and Earth
orientation inputs. For propagation record the TLE, gravity convention, offsets
and output frame. Compare the same scientific output with the same inputs.
When a backend changes model, time coverage or accuracy, state the trade-off
alongside its latency.

For warmed calls, prepare stable fixtures, contexts and reusable objects before
`b.iter`. For end-to-end searches, include any cache construction or allocation
the API performs on each call. Criterion's `iter` consumes the returned value
through a compiler barrier and includes its destruction. Place input barriers on
values whose compile-time visibility could change the work and return the result
from the closure. Verify expensive paths in a profile or optimized assembly
when the observed cost looks implausible.

Use `iter_batched` or `iter_batched_ref` when each timed operation needs fresh
untimed setup. Document batch size, retained values and cache effects. For a
genuinely cold load, define process and OS page-cache state separately and
measure the load as its own operation. A parser given preloaded bytes measures
CPU and allocation for parsing. A warm repeated fixed-epoch query measures
latency with favorable record locality. A realistic streaming workload needs a
bounded, deterministic epoch trace or an explicitly defined batch. Report
batch latency and distinguish it from throughput computed per item.

## Capture the environment

Start from the repository root. Record the source revision and dirty diff,
toolchain, exact command and UTC measurement date. Keep the lockfile and resolved
feature graph with the results. These are example capture commands for Linux.

```sh
git rev-parse HEAD
git status --short
git diff HEAD
rustc -vV
cargo -V
cargo tree --locked -p siderust --features bench-internals -e features
uname -a
lscpu
date -u
```

Record effective target triple, optimization, LTO, codegen units, allocator,
linker, CPU features (including any `target-cpu=native`), Cargo configuration
and environment overrides such as `RUSTFLAGS` and `CARGO_PROFILE_BENCH_*`.
For external files record origin/version, `sha256sum`, coverage, and kernel
load order where applicable. Include active time/EOP data identity and
observed or predicted status when the workload uses those tables.

Use the same machine, power mode and affinity for the baseline and candidate.
Record governor, boost, SMT/core placement and relevant cache and memory
configuration. Minimize competing CPU and I/O activity and check sustained
temperature and frequency. Containers, virtual machines and shared runners
can add run-to-run variance even when one run looks stable.

## Measure and compare

Use the target's configured Criterion settings for an initial run. In the
current suites, solar, lunar and stellar searches request 20 samples and 5 s
of measurement. Cross-body altitude requests 20 and 10 s. Converters request
1,000 and 3 s with 1 s warm-up. SGP4 requests 30 and 2 s with 1 s warm-up.
The ephemeris comparison requests 100 and 5 s. Other settings inherit
Criterion defaults (including 3 s warm-up, 95% confidence and 5% significance).
Measurement time is a target. Long operations may exceed
it to satisfy the sample count. For a consistently expensive operation,
consider flat sampling so later samples do not grow with the sample index.
Record any overrides and apply them to both sides before collecting the
comparison.

```sh
cargo bench --locked -p siderust --bench sgp4 -- \
  --exact sgp4/leo/propagate --save-baseline before
# Run the candidate with the same harness, inputs and build configuration.
cargo bench --locked -p siderust --bench sgp4 -- \
  --exact sgp4/leo/propagate --baseline before
```

Use a fresh baseline name per experiment. Preserve the named baseline under
the Criterion output directory so both builds can use it. Use a common
absolute `CRITERION_HOME` to share results between separate build directories. Across revisions,
check the harness diff and resolved dependencies first. If the harness changed,
use the same harness with both revisions when their interfaces permit it.
Record the source of that harness and any patches. Treat workload or model
changes as a fresh benchmark series even when case IDs survive. If the lockfile
changed deliberately, report whether the experiment measures the whole revision
or holds dependencies fixed.

Criterion warms up the case, times samples containing multiple iterations and
uses resampling to estimate uncertainty within the run. Keep the estimator,
time units, confidence level and interval in the report. The interval does not
capture variation between machines or independent runs. Repeat measurements
with alternating baseline and candidate order when practical, especially for
small effects. Statistical significance describes evidence of a difference
under the sampling assumptions. The practical effect depends on the absolute
cost, workload frequency and acceptable budget. Investigate inconsistent
replicates and unexplained outliers before calling a regression.

## Validate and report

Run the relevant [scientific tests](coverage.md) alongside performance
measurements. Where a comparison has a correctness preflight, run it before
sampling. Record fixture and tolerance, compare event counts and directions
as well as timestamps, and include coverage or failure behavior if it informs
backend choice. Finite results and successful benchmark execution establish
only that the workload completed.

Keep a short report with the question, exact workload, configuration, reference,
time estimates and uncertainty, interpretation, limitations and reproduction
commands. Archive small logs and Criterion baseline directories with the
experiment when needed. Criterion already stores machine-readable estimates
and samples. Add a separate schema only when those artifacts cannot support
an established consumer. Keep generated trees out of the source checkout.

CI type-checks the benchmark targets and checks relative documentation links.
Use `cargo bench --no-run` locally to verify optimized codegen and linking,
and Criterion `--test` to verify representative executable cases. With
compact-only features, select individual `--bench` targets. The plural
`--benches` also compiles the library unit-test harness, which currently
assumes default library features. Use a
[profile](profiling.md) to investigate reproducible changes before collecting
fresh uninstrumented timings.
