## 1. What Flamegraphs Are (and What They Are Not)

A flamegraph is a **statistical CPU profile visualization**:

* The **width** of a frame represents how much CPU time was spent in that function (inclusive).
* The **vertical axis** represents the call stack depth.
* The **top frames** are where the CPU is actually executing.

Flamegraphs answer:

* *Where is the CPU time going?*
* *Which functions dominate runtime?*

They do **not**:

* Measure wall‑clock latency precisely
* Explain memory allocation behavior (use heap profilers for that)
* Replace Criterion’s statistical benchmarking

Criterion + flamegraphs is the correct combination: **Criterion detects regressions, flamegraphs explain them**.

---

## 2. System Requirements (Linux)

### 2.1 Required tools

```bash
cargo install flamegraph
```

You also need:

* `perf` (Linux performance tools)
* Kernel support for performance counters

On Debian/Ubuntu:

```bash
sudo apt install linux-tools-common linux-tools-$(uname -r)
```

On Arch:

```bash
sudo pacman -S perf
```

---

## 3. Kernel Configuration (Critical)

Most modern distributions restrict access to `perf` by default. If you see errors mentioning `perf_event_paranoid`, you must relax the kernel settings.

### 3.1 Temporary (recommended for development)

```bash
sudo sysctl -w kernel.perf_event_paranoid=1
sudo sysctl -w kernel.kptr_restrict=0
```

This change lasts until reboot.

### 3.2 Permanent (optional)

Create:

```bash
sudo nano /etc/sysctl.d/99-perf.conf
```

Add:

```conf
kernel.perf_event_paranoid = 1
kernel.kptr_restrict = 0
```

Apply:

```bash
sudo sysctl --system
```

### 3.3 Alternative: capabilities instead of sysctl

If you cannot change system settings globally:

```bash
sudo setcap cap_perfmon,cap_sys_ptrace=eip $(which perf)
```

---

## 4. Rust Build Configuration for Profiling

Flamegraphs are only useful if symbols are available.

In `Cargo.toml`:

```toml
[profile.release]
debug = true
```

Notes:

* This **does not disable optimizations**.
* It only adds DWARF symbols so stacks and function names are visible.

---

## 5. Profiling Criterion Benchmarks

### 5.1 Basic command

```bash
cargo flamegraph --bench solar_altitude
```

This:

* Builds the benchmark in `release`
* Runs it under `perf record`
* Produces `flamegraph.svg`

Open it with any browser.

---

### 5.2 Profiling a specific Criterion benchmark

Criterion supports filtering benchmarks by name. Arguments after `--` are passed to the benchmark binary.

Example:

```bash
cargo flamegraph --bench solar_altitude -- --bench find_night_periods_3
```

Use this to avoid profiling unrelated benchmarks.

---

## 6. Reducing Noise from Criterion

Criterion performs:

* Warm‑up iterations
* Multiple samples
* Statistical aggregation

This means flamegraphs may include **more runtime than strictly necessary**.

### 6.1 Recommended approach

For deep profiling:

* Create a **dedicated benchmark** for the code path you want to analyze
* Reduce noise:

```rust
Criterion::default()
    .sample_size(10)
    .warm_up_time(Duration::from_secs(1))
```

### 6.2 Preventing dead‑code elimination

Always protect inputs and outputs:

```rust
use criterion::black_box;

black_box(your_function(black_box(input)));
```

If you do not do this, the optimizer may remove the work entirely, producing misleading flamegraphs.

---

## 7. Interpreting Flamegraphs Correctly

### 7.1 What to look for

* Very wide frames near the top → **primary optimization targets**
* Repeated deep stacks → abstraction or iterator overhead
* Unexpected functions → hidden costs (allocations, formatting, bounds checks)

### 7.2 What *not* to optimize prematurely

* Thin frames (<5–10%)
* Code that is already memory‑bound
* Setup/teardown outside the hot path

Optimize the **widest stack first**.

---

## 8. Common Pitfalls

* ❌ Profiling in `debug` mode
* ❌ Forgetting `debug = true` in release
* ❌ Interpreting flame height as time (only width matters)
* ❌ Optimizing without re‑benchmarking in Criterion

Always:

1. Measure with Criterion
2. Explain with flamegraph
3. Optimize
4. Re‑measure

---

## 9. When `perf` Is Not Available

If system policy prevents `perf`, use `samply`:

```bash
cargo install samply
samply record cargo bench --bench solar_altitude -- --bench find_night_periods_3
```

`samply` provides an interactive profiler UI and is often easier to use on locked‑down systems.

---

## 10. Recommended Workflow for Siderust

1. Detect regression with Criterion
2. Run `cargo flamegraph` on the affected benchmark
3. Identify top 1–2 hot functions
4. Optimize algorithmically first
5. Re‑run Criterion
6. Commit with benchmark evidence

This keeps performance work **intentional, explainable, and reproducible**.

---

## 11. Appendix: Useful Commands

```bash
# Check kernel perf settings
sysctl kernel.perf_event_paranoid kernel.kptr_restrict

# Verify perf installation
perf --version

# Clean old profiling artifacts
rm -f flamegraph.svg perf.data
```
