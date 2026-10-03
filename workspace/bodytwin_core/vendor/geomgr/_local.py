"Adapters for optional workstation geometry and export dependencies."
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[3]
RES = WS / 'results'
BT_SRC = WS / 'references' / 'current_bodytwin' / 'src'
MSK_PY = Path('source_repository/.venv-msk/bin/python')            # OpenSim 4.6
FIELD_PY = Path('~/projects/3fold_staging/.venv-bodytwin/bin/python')  # warp (CPU) for the field engine
STAGING = Path('~/projects/3fold_staging')

sys.dont_write_bytecode = True


def _path(p):
    p = str(p)
    if p not in sys.path:
        sys.path.insert(0, p)


def bodytwin_geometry_certificate(V, F):
    _path(BT_SRC)
    from bodytwin.geometry.geometry_validity_v1 import geometry_certificate
    return geometry_certificate(V, F, evaluate=True)


_CONTRACT = None


def contract():
    """BT-C1b contract module loaded under a private name (read-only)."""
    global _CONTRACT
    if _CONTRACT is None:
        f = RES / 'BT-C1b' / 'compose' / 'contract.py'
        spec = importlib.util.spec_from_file_location('n7a_bt_c1b_contract', f)
        mod = importlib.util.module_from_spec(spec)
        sys.modules['n7a_bt_c1b_contract'] = mod
        spec.loader.exec_module(mod)
        _CONTRACT = mod
    return _CONTRACT


