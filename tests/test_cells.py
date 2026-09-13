"""Every cell runs standalone on CPU and exits 0, writing only under BODYTWIN_OUT.

Cells that read another cell's result JSON run after their producers, all into one shared
BODYTWIN_OUT; producers that are not part of this repository are replaced by the synthetic
inputs under examples/synthetic_inputs/ (copied in before the first run).
"""
import json
import os
import shutil
import subprocess
import sys
import time

import pytest

from conftest import EXAMPLES, REPO, cell_scripts

SCRIPTS = {p.stem: p for p in cell_scripts()}
DEPS = json.loads((REPO / "tests" / "cell_dependencies.json").read_text())
SYNTHETIC = EXAMPLES / "synthetic_inputs"


def run_order():
    order, seen = [], set()

    def visit(name):
        if name in seen:
            return
        seen.add(name)
        for dep in DEPS.get(name, []):
            if dep in SCRIPTS:
                visit(dep)
        order.append(name)

    for name in sorted(SCRIPTS):
        visit(name)
    return order


ORDER = run_order()


DEFAULT_TIMEOUT = 120
# Cells whose own selftest is longer than the default cap; they have no size/iteration flag to
# reduce, so they get a wider cap instead of a weakened run.
TIMEOUTS = {
    "cellcycle_restriction_switch_biomd265_rebuild": 400,
    "fingertip_tactile": 400,
    "gi_motility_slow_waves": 400,
}


class CellRunner:
    """Runs a cell (and, first, its producers) once and caches the result.

    Lazy so that a filtered run (pytest -k) executes only the selected cells and the producers
    they read, not the whole tree.
    """

    def __init__(self, out):
        self.out = out
        self.env = dict(os.environ)
        self.env["BODYTWIN_OUT"] = str(out)
        self.cache = {}
        self.running = set()

    def __getitem__(self, name):
        if name in self.cache:
            return self.cache[name]
        if name in self.running:
            raise AssertionError(f"cyclic dependency at {name}")
        self.running.add(name)
        for dep in DEPS.get(name, []):
            if dep in SCRIPTS:
                self[dep]
        script = SCRIPTS[name]
        t0 = time.time()
        try:
            proc = subprocess.run([sys.executable, script.name], cwd=script.parent, env=self.env,
                                  capture_output=True, text=True,
                                  timeout=TIMEOUTS.get(name, DEFAULT_TIMEOUT))
            result = (proc.returncode, proc.stdout, proc.stderr, time.time() - t0)
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout or ""
            if isinstance(stdout, bytes):
                stdout = stdout.decode(errors="replace")
            result = (None, stdout, f"timeout after {TIMEOUTS.get(name, DEFAULT_TIMEOUT)} s",
                      time.time() - t0)
        self.running.discard(name)
        self.cache[name] = result
        return result


@pytest.fixture(scope="session")
def cell_runs(tmp_path_factory):
    out = tmp_path_factory.mktemp("bodytwin_out")
    if SYNTHETIC.is_dir():
        for producer in SYNTHETIC.iterdir():
            if producer.is_dir():
                shutil.copytree(producer, out / producer.name, dirs_exist_ok=True)
    return CellRunner(out)


@pytest.mark.parametrize("name", ORDER, ids=[SCRIPTS[n].parent.name + "/" + n for n in ORDER])
def test_cell_runs(name, cell_runs):
    rc, stdout, stderr, _ = cell_runs[name]
    assert rc == 0, stdout[-2000:] + stderr[-2000:]
    if name != "ventricular_ap_core_tt04style":  # a shared solver library, not a cell
        assert stdout.strip(), "a cell must print its own numbers"


def test_every_cell_writes_only_under_out(cell_runs):
    stray = [str(p.relative_to(REPO)) for p in (REPO / "src").rglob("*.json")]
    assert stray == [], stray


def test_every_subsystem_has_at_least_one_cell():
    subsystems = {p.parent.name for p in SCRIPTS.values()}
    assert len(subsystems) >= 15


def test_synthetic_inputs_are_marked():
    if not SYNTHETIC.is_dir():
        pytest.skip("no synthetic inputs")
    for path in SYNTHETIC.rglob("*.json"):
        data = json.loads(path.read_text())
        assert data.get("synthetic") is True, path
