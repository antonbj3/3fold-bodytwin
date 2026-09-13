import os
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
CELLS = REPO / "src" / "bodytwin" / "cells"
FRAMEWORK = REPO / "src" / "bodytwin" / "framework"
EXAMPLES = REPO / "examples"

sys.path.insert(0, str(FRAMEWORK))


def cell_scripts():
    return sorted(p for p in CELLS.rglob("*.py") if p.name != "__init__.py")


def framework_env(root, graph=None):
    env = dict(os.environ)
    env["BODYTWIN_ROOT"] = str(root)
    env["BODYTWIN_GRAPH"] = str(graph or pathlib.Path(root) / "ANCHOR_GRAPH.json")
    env["BODYTWIN_INVENTORY"] = str(pathlib.Path(root) / "CELL_INVENTORY.md")
    return env
