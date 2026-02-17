# Architecture Docs

This folder is intentionally **not** a dumping ground for reports.
It is meant to capture the project’s **architectural decisions** and the
**promises/guarantees** that users can rely on across releases.

If you are looking for benchmarking/profiling write-ups, see `benches/reports/README.md`.

## Documents

- `doc/architecture_decisions.md`: key architectural decisions (what/why/impact).
- `doc/promises.md`: compatibility/behavior promises (what we guarantee).
- `doc/datasets.md`: how embedded datasets and optional JPL backends work.

## Related docs elsewhere

- `README.md`: user-facing overview and installation.
- `examples/README.md`: runnable tour of features.
- `src/coordinates/README.md`: coordinate architecture and invariants (module-level).
- `benches/README.md`: benchmark entrypoints and how to run them.
