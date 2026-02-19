# dev-deps

This folder contains **git submodules** (`affn`, `cheby`, `qtty`, `tempoch`) that are kept in this repository
**only for development and CI convenience** (local hacking, debugging, and testing across the Siderust crate
ecosystem).

## Publishing rule

The released/published `siderust` crate must **always** depend on the **published crates from crates.io**
(as declared in `Cargo.toml`) and must **not** switch to `path = "..."`
dependencies pointing at `dev-deps/*`.

## Optional local overrides (do not commit)

If you want to test `siderust` against local changes in these submodules, use a local patch override in your
working copy (and do not commit it), e.g.:

```toml
# Cargo.toml (local only)
[patch.crates-io]
affn = { path = "dev-deps/affn" }
cheby = { path = "dev-deps/cheby" }
qtty = { path = "dev-deps/qtty" }
tempoch = { path = "dev-deps/tempoch" }
```
