"""Verify the frozen research contract and executable release before replay."""
import hashlib, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def check():
    lock = ROOT / 'DELIVERY_CODE_LOCK.json'
    if sha(lock) != lock.with_suffix('.sha256').read_text().strip():
        raise ValueError('Release lock corrupt')
    for (rel, h) in json.loads(lock.read_text())['files'].items():
        if sha(ROOT / rel) != h:
            raise ValueError('Release code drift: ' + rel)
    for p in list(ROOT.glob('PREREG_*.json')) + list(ROOT.glob('FROZEN_PREDICTIONS*.json')):
        if sha(p) != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('Research freeze drift: ' + p.name)
    print('Research freezes and executable code lock verified')
if __name__ == '__main__':
    check()
