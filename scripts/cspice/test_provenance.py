"""Offline checks for invalid generator configuration and kernel provenance."""

import hashlib
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from provenance import check_kernel, check_toolkit, provenance


class ProvenanceTests(unittest.TestCase):
    def test_missing_and_corrupt_kernel(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "de440s.bsp"
            with self.assertRaisesRegex(ValueError, "missing kernel"):
                check_kernel(path, path.name)
            path.write_bytes(b"invalid")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                check_kernel(path, path.name)

    def test_matching_kernel(self):
        content = b"kernel contents"
        digest = hashlib.sha256(content).hexdigest()
        with TemporaryDirectory() as directory:
            path = Path(directory) / "example.bsp"
            path.write_bytes(content)
            with patch(
                "provenance.provenance",
                return_value={"source_sha256": {path.name: digest}},
            ):
                self.assertEqual(check_kernel(path, path.name), digest)

    def test_toolkit_mismatch(self):
        expected = provenance()
        version = expected["generator"].removeprefix("SpiceyPy ")
        check_toolkit(
            SimpleNamespace(__version__=version, tkvrsn=lambda _: expected["toolkit"])
        )
        with self.assertRaisesRegex(ValueError, "expected"):
            check_toolkit(
                SimpleNamespace(__version__=version, tkvrsn=lambda _: "other")
            )

    def test_invalid_output_mode_exits_before_generation(self):
        script = Path(__file__).with_name("cspice_lunar_de440_references.py")
        with TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, script, "--verify", "--output-dir", directory],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("not allowed with argument", result.stderr)
