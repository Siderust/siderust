# Maintain the documentation

Use the root README for project orientation and navigation. Use the
manual for tutorials, tasks, physical explanations, and cross-cutting
contracts. Put signatures, input/output semantics, errors, panics, and
small runnable API examples in Rustdoc. Keep complete workflows in
`examples/`, reference comparisons near their fixtures, and performance
measurements in [benchmark reports](../../benches/reports/README.md).
Use the changelog for release history and issues or pull requests for
temporary implementation coordination.

When documenting a scientific API, state the target, center, frame,
units, time scale, corrections, and domain where they affect the result.
Cite the model or convention and identify the independent reference used
for numerical comparisons. Record the software revision, input provenance,
test conditions, and tolerance with measured results.

At API boundaries, prefer typed quantities and epochs. Reuse `qtty`
types for general physical quantities and define astronomy-specific concepts
in Siderust. Keep conversion from external scalar data near its boundary.
Prefer marker traits for compile-time model selection. Isolate runtime
dispatch at interfaces that load or select models at runtime.

Run these checks from the repository root for changes to links, API docs,
or examples. Repeat example and test checks with the features used in
the documented command. The [CI workflow](../../.github/workflows/github-ci.yml)
contains the supported build matrix.

```sh
cargo run --locked -p siderust-release --bin check-doc-links
RUSTDOCFLAGS="-D rustdoc::broken_intra_doc_links -D rustdoc::invalid_html_tags" cargo doc -p siderust --all-features --no-deps
cargo test -p siderust --doc --all-features
cargo check -p siderust --examples --all-features
```

Run tutorials and runnable examples after editing them. Inspect Markdown
headings, tables, and code fences in the rendered page when adding a new
pattern. Update inbound links and anchors when moving content. Leave
historical filenames in changelog prose when they describe the release
at the time.
