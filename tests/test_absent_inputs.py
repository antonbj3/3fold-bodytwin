"""Absent-input handling for three cells (a fail-fast fix from the 2026-09 review).

A cell whose REQUIRED upstream result is missing must stop with exit 1 and one
``FAIL: required input missing`` line instead of a raw NameError/UnboundLocalError.
pulmonary_gas_exchange's respiratory input is OPTIONAL: with it absent the cell must
not raise a TypeError while formatting the missing value.
"""
import os
import shutil
import subprocess
import sys

from conftest import CELLS, EXAMPLES

RAW_EXCEPTIONS = ("NameError", "UnboundLocalError", "TypeError", "Traceback")


def _run(rel, out, timeout=300):
    script = CELLS / rel
    env = dict(os.environ, BODYTWIN_OUT=str(out))
    env.pop("BODYTWIN_I1_REGISTRY", None)
    env.pop("BODYTWIN_EVIDENCE_MODE", None)
    return subprocess.run([sys.executable, script.name], cwd=script.parent, env=env,
                          capture_output=True, text=True, timeout=timeout)


def _assert_required_missing(proc):
    text = proc.stdout + proc.stderr
    assert proc.returncode == 1, text[-2000:]
    assert "FAIL: required input missing" in proc.stdout, text[-2000:]
    for exc in RAW_EXCEPTIONS:
        assert exc not in text, text[-2000:]


def test_muscle_energetics_required_inputs_missing(tmp_path):
    _assert_required_missing(_run("energy/muscle_energetics.py", tmp_path))


def test_thyroid_metabolic_axis_required_input_missing(tmp_path):
    _assert_required_missing(_run("endocrine/thyroid_metabolic_axis.py", tmp_path))


def test_pulmonary_gas_exchange_required_input_missing(tmp_path):
    proc = _run("respiratory/pulmonary_gas_exchange.py", tmp_path)
    text = proc.stdout + proc.stderr
    assert proc.returncode == 1, text[-2000:]
    for exc in RAW_EXCEPTIONS:
        assert exc not in text, text[-2000:]


def test_pulmonary_gas_exchange_optional_respiratory_input_absent(tmp_path):
    # The required cardiac_output result is produced from the committed synthetic inputs;
    # the optional respiratory result is deliberately not produced.
    for producer in (EXAMPLES / "synthetic_inputs").iterdir():
        if producer.is_dir():
            shutil.copytree(producer, tmp_path / producer.name, dirs_exist_ok=True)
    cardiac = _run("cardiovascular/cardiac_output.py", tmp_path)
    assert (tmp_path / "cardiac_output" / "cardiac_output_results.json").exists(), \
        (cardiac.stdout + cardiac.stderr)[-2000:]
    assert not (tmp_path / "respiratory" / "respiratory_results.json").exists()
    proc = _run("respiratory/pulmonary_gas_exchange.py", tmp_path)
    text = proc.stdout + proc.stderr
    for exc in RAW_EXCEPTIONS:
        assert exc not in text, text[-2000:]
    assert "rest=n/a mL/min" in proc.stdout, text[-2000:]
