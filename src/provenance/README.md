# `provenance`

## Scientific scope

Reproducible astronomical software must preserve the origin of every dataset it
uses. Provenance is the metadata that answers where a table, spectrum, kernel,
or derived product came from, which release or citation it corresponds to, and
which exact bytes were consumed. This module provides the metadata vocabulary
used in `siderust` to make those questions explicit.

The intent is not to model the full workflow graph of an observatory archive.
Instead, the focus is the pragmatic provenance needed by a numerical library:
source identity, version, retrieval or generation timestamp, integrity hash, and
notes about assumptions or post-processing.

## Technical scope

What you will find here:

- `Provenance`: the main record type.
- `DataSource`: origin classifier for literature, bundled files, external
  resources, and computed products.
- Builder-style constructors and enrichers such as `bundled_file`, `cited`,
  `computed`, `with_version`, and `with_notes`.
- `checksum.rs`: compile-time SHA-256 support and checksum-assertion helpers.

This module is also re-exported by feature-gated data containers such as
`spectra` and `tables`, so the same provenance vocabulary is shared across the
crate.

## References

- International Virtual Observatory Alliance. Provenance Data Model, version
  1.0.
  <https://www.ivoa.net/documents/ProvenanceDM/>
- IERS Conventions (2010), IERS Technical Note 36, Chapter 1.
  <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36>
