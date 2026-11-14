"""Validate the reference environment against the published DE440 provenance."""

import hashlib
import json
from pathlib import Path

MANIFEST = (
    Path(__file__).resolve().parents[2] / "tests/fixtures/lunar_de440/provenance.json"
)


def provenance():
    return json.loads(MANIFEST.read_text(encoding="ascii"))


def check_kernel(path, name):
    expected = provenance()["source_sha256"][name]
    if not path.is_file():
        raise ValueError(f"missing kernel: {path}")
    with path.open("rb") as kernel:
        actual = hashlib.file_digest(kernel, "sha256").hexdigest()
    if actual != expected:
        raise ValueError(
            f"SHA-256 mismatch for {path}: expected {expected}, got {actual}"
        )
    return actual


def check_kernels(directory):
    return {
        name: check_kernel(directory / name, name)
        for name in provenance()["source_sha256"]
    }


def check_toolkit(spice):
    expected = provenance()
    generator = f"SpiceyPy {spice.__version__}"
    toolkit = spice.tkvrsn("TOOLKIT")
    if (generator, toolkit) != (expected["generator"], expected["toolkit"]):
        raise ValueError(
            f"expected {expected['generator']} / {expected['toolkit']}, got {generator} / {toolkit}"
        )
