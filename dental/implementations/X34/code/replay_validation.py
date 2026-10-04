"""Replay frozen R2B operations; serialize deliberately injected NaN as a labelled null.

The original frozen harness correctly rejected NaN, then tried to write it into
strict JSON in a test DETAIL. Keep that code/failure; no numerical rule is changed.
"""
import math
from pathlib import Path
import json
import validate_r2b as V
from initialize import H, sha

def clean(x):
    if isinstance(x, float) and (not math.isfinite(x)):
        return {'injected_nonfinite': repr(x)}
    if isinstance(x, dict):
        return {k: clean(v) for (k, v) in x.items()}
    if isinstance(x, list):
        return [clean(v) for v in x]
    return x

def safe_dump(p, x):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def main():
    lock = json.loads((H / 'REPLAY_LOCK.json').read_text())
    for (p, h) in lock['code'].items():
        if sha(H / p) != h:
            raise ValueError('Replay code drift: ' + p)
    V.dump = safe_dump
    return V.main()
if __name__ == '__main__':
    raise SystemExit(main())
