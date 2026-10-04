from pathlib import Path
import hashlib, json
P = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json', 'DECOMPOSITION_R1.json', 'DECOMPOSITION_R2.json', 'DECOMPOSITION_R3.json', 'DECOMPOSITION_R4.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R3.json']:
    p = P / name
    assert p.is_file(), name + ' missing'
    assert sha(p) == p.with_suffix('.json.sha256').read_text().strip(), name + ' hash mismatch'
for row in json.loads((P / 'REVIEW_LOCAL_LOCK.json').read_text())['files']:
    q = (P / row['path']).resolve()
    assert q.is_relative_to(P), 'Nonlocal execution input'
    assert sha(q) == row['sha256'], 'Local code/input drift: ' + row['path']
print('Frozen and actually executed local code/input hashes verified')
