"""All numerical controls run in code/validate.py and its saved per-run raw JSON.

They compare external published geometry, exact unchanged GenCAD witnesses,
injected bad inputs and transport corruption. This file exposes the replay to
standard unittest without duplicating implementation-shaped assertions.
"""
import subprocess
import unittest
from pathlib import Path

class FrozenCapability(unittest.TestCase):

    def test_frozen_numerical_and_external_controls(self):
        root = Path(__file__).resolve().parents[1]
        self.assertEqual(subprocess.run([str(root / 'run_all.sh')], cwd=root).returncode, 0)
