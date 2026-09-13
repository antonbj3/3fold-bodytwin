"""The coupled metabolic chain runs standalone on CPU, exits 0, and writes only under BODYTWIN_OUT."""
import json
import os
import subprocess
import sys

from conftest import REPO

CHAIN = REPO / "src" / "bodytwin" / "chains" / "metabolic_chain_v1.py"


def test_metabolic_chain_v1(tmp_path):
    env = dict(os.environ, BODYTWIN_OUT=str(tmp_path))
    proc = subprocess.run([sys.executable, str(CHAIN)], env=env, cwd=str(REPO),
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    results = json.loads((tmp_path / "metabolic_chain_v1" /
                          "metabolic_chain_v1_results.json").read_text())
    assert results["overall_pass"] is True
    assert all(results["gates"].values())
    assert results["n_draws"] >= 64
    assert results["propagation"]["full_chain"]["half_width_l_min"] <= \
        results["propagation"]["linear_sum_of_stage_half_widths_l_min"]
    assert results["propagation"]["nonlinear_stage_flagged"] is None
