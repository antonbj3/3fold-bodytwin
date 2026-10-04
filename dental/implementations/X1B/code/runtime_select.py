"""Stdlib-only launcher: use the explicitly pinned installed scientific runtime."""
from pathlib import Path
import json
print(json.loads((Path(__file__).resolve().parents[1] / 'RUNTIME_LOCK.json').read_text())['interpreter_locator'])
