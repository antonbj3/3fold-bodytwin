"""Freeze once, verify on rerun; never silently replace a prediction artifact."""
import hashlib, json
from pathlib import Path

def frozen_write(path, payload):
    path = Path(path)
    if path.exists():
        old = json.loads(path.read_text())
        a = dict(old)
        b = dict(payload)
        for x in [a, b]:
            x.pop('created_at_utc', None)
            x.pop('prediction_payload_sha256', None)
        if a != b:
            raise ValueError(f'frozen prediction changed: {path}; archive explicitly and version before continuing')
        return
    stable = {k: v for (k, v) in payload.items() if k != 'created_at_utc'}
    payload['prediction_payload_sha256'] = hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
