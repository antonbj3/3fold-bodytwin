"""Regression for the chain boundary invariants (from the 2026-09 review).

A wrong-unit swap in stage 2 (O2 reported in L/min in a mL/min slot, names and keys
unchanged) moves the observable while g01-g06 still pass; g07 is the gate that must
catch it (O2/glucose about 0.0056 mol/mol instead of about 5.6).
"""
import json
import os
import subprocess
import sys
import textwrap

from conftest import REPO

CHAIN = REPO / "src" / "bodytwin" / "chains" / "metabolic_chain_v1.py"

DRIVER = textwrap.dedent("""
    import contextlib, importlib.util, io, json, sys
    spec = importlib.util.spec_from_file_location("chain_under_test", sys.argv[1])
    m = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(m)
    swap = sys.argv[2] == "swap"
    if swap:
        correct = m.stage2_atp_and_o2
        def swapped(disposal_mg_min, po_nadh):
            out = correct(disposal_mg_min, po_nadh)
            out["o2_required_ml_min"] = out["o2_required_ml_min"] / 1000.0
            return out
        m.stage2_atp_and_o2 = swapped
    with contextlib.redirect_stdout(io.StringIO()):
        rc = m.main()
    results = json.loads(m.OUT_JSON.read_text())
    print(json.dumps({"rc": rc, "gates": results["gates"],
                      "overall_pass": results["overall_pass"],
                      "boundary_invariants": results["boundary_invariants"],
                      "q_l_min": results["nominal"]["final_cardiac_output_l_min"]}))
""")


def _run(tmp_path, variant):
    env = dict(os.environ, BODYTWIN_OUT=str(tmp_path / variant))
    proc = subprocess.run([sys.executable, "-c", DRIVER, str(CHAIN), variant], env=env,
                          cwd=str(REPO), capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return json.loads(proc.stdout.strip().splitlines()[-1])


def test_correct_units_pass_all_gates_including_g07(tmp_path):
    out = _run(tmp_path, "correct")
    assert out["rc"] == 0 and out["overall_pass"] is True
    assert out["gates"]["g07_boundary_invariants_hold"] is True
    ratio = out["boundary_invariants"]["o2_per_glucose_mol_per_mol"]
    assert 4.0 <= ratio <= 8.0


def test_wrong_unit_swap_is_caught_only_by_g07(tmp_path):
    good = _run(tmp_path, "correct")
    bad = _run(tmp_path, "swap")
    gates = bad["gates"]
    for name, value in gates.items():
        if name != "g07_boundary_invariants_hold":
            assert value is True, (name, gates)
    assert gates["g07_boundary_invariants_hold"] is False
    assert bad["overall_pass"] is False and bad["rc"] == 1
    assert bad["boundary_invariants"]["o2_per_glucose_in_band"] is False
    assert bad["boundary_invariants"]["o2_per_glucose_mol_per_mol"] < 0.01
    assert bad["q_l_min"] < good["q_l_min"]


PARAMS_DRIVER = textwrap.dedent("""
    import contextlib, importlib.util, json, sys, io
    spec = importlib.util.spec_from_file_location("chain_params_under_test", sys.argv[1])
    m = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(m)
    doubled = dict(m.SYNTHETIC_REFERENCE,
                   meal_glucose_load_mg=2.0 * m.SYNTHETIC_REFERENCE["meal_glucose_load_mg"])
    run = m.run_chain(doubled)
    print(json.dumps({"with_params": m.boundary_invariants(run, doubled),
                      "default_params": m.boundary_invariants(run)}))
""")


def test_boundary_invariants_use_the_parameters_of_the_run(tmp_path):
    env = dict(os.environ, BODYTWIN_OUT=str(tmp_path))
    proc = subprocess.run([sys.executable, "-c", PARAMS_DRIVER, str(CHAIN)], env=env,
                          cwd=str(REPO), capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = json.loads(proc.stdout.strip().splitlines()[-1])
    # a doubled meal closes its mass balance only with its own parameters
    assert out["with_params"]["glucose_mass_balance"] is True
    assert out["with_params"]["fick_closure"] is True
    assert out["default_params"]["glucose_mass_balance"] is False
