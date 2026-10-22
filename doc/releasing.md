# Releasing Siderust

Maintainer checklist for publishing a versioned release (for example `v0.15.0`)
to crates.io and GitHub.

## Prerequisites

1. Merge the feature / release-prep PRs into `main` (including any stacked PRs
   that must land first).
2. Confirm CI is green on `main`.
3. Confirm public dependencies resolve from crates.io (no temporary Git/path
   patches in publishable packages).
4. Confirm `CHANGELOG.md` has a nonempty `## [X.Y.Z]` section (date optional
   until the tag is cut) and that `siderust` / `siderust-ffi` report `X.Y.Z`
   (FFI remains `publish = false` unless intentionally opened).

Local preflight (no registry token required):

```bash
cargo test -p siderust-release
cargo build -p siderust-release
target/debug/siderust-release validate
bash scripts/publish-changed.sh --dry-run
```

## Cut the release

1. On the commit that should ship, create and push the annotated or lightweight
   version tag:

   ```bash
   git tag v0.15.0
   git push origin v0.15.0
   ```

2. The **Publish to crates.io** workflow runs on that tag:
   - validates tag ↔ Cargo ↔ changelog consistency;
   - publishes publishable crates via `scripts/publish-changed.sh`;
   - records a GitHub deployment against the `crates.io` environment;
   - creates (or updates) a GitHub Release titled `Siderust vX.Y.Z` with notes
     extracted from `CHANGELOG.md`.

3. Verify:
   - [crates.io/crates/siderust](https://crates.io/crates/siderust) shows the
     new version;
   - the repository **Deployments** entry for `crates.io` succeeded;
   - the GitHub Release notes match the changelog section.

Do **not** create the GitHub Release manually on ordinary commits or PRs; only
the tag-triggered workflow should publish notes after a successful crate
upload.

## Recovery

| Situation | Action |
|-----------|--------|
| Validation fails before publish | Fix metadata/changelog on a follow-up commit; delete the bad tag if it was pushed early, retag the fixed commit, push again. |
| Publish fails mid-job | Fix the cause; re-run the failed workflow. Already-uploaded versions are skipped via the crates.io version API, so reruns are idempotent. |
| Publish succeeded, GitHub Release failed | Re-run the `github-release` job (or the whole workflow). Existing releases are updated in place; crates are not republished. |
| Wrong notes on an existing release | Fix `CHANGELOG.md` on `main` if needed, then re-run the release job for the same tag, or `gh release edit`. |
| Accidental publish of a bad crate version | Yank on crates.io if appropriate; ship a patch release. Do not reuse the same version number. |

## FFI crate

`siderust-ffi` stays `publish = false` by default. To publish it deliberately:

1. Set `publish = true` in `siderust-ffi/Cargo.toml`.
2. Run `bash scripts/publish-changed.sh --confirm-ffi` with a registry token
   (prefer a protected environment).
3. Revert `publish` to `false` unless the FFI crate should remain publishable.
