"""Explicit verifier-only revision; original freezes/code/failure are archived."""
import sys, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from seal import write_anchor
A = ROOT / 'revisions/RELEASE_1'
old = read(A / 'BENCHMARK_LOCK.json')
previous = read(A / 'FROZEN_PREDICTIONS.json')
assert sha(ROOT / 'BENCHMARK_LOCK.json') == sha(A / 'BENCHMARK_LOCK.json')
allowed = {'code/scorer.py', 'code/controls.py', 'code/legacy_v1/io.py', 'SOURCE_REUSE.json'}
changed = []
for (name, entry) in old['files'].items():
    actual = sha(ROOT / name)
    if actual != entry['sha256']:
        if name not in allowed:
            raise RuntimeError('unexpected input change ' + name)
        changed.append(name)
assert set(changed) == allowed
for (rel, info) in previous['files'].items():
    assert sha(PAYLOAD / 'predictions' / rel) == info['sha256']
new = read(A / 'BENCHMARK_LOCK.json')
new.update(frozen_utc=now(), revision='R8_VERIFIER_CORRECTION', parent_benchmark_sha256=sha(A / 'BENCHMARK_LOCK.json'), revision_scope='Only workspace-independent legacy IO and fail-closed metric arithmetic; no source, task, thresholds, generator, runtime or prediction array changes')
for name in changed + ['PREREG_R8.json']:
    p = ROOT / name
    new['files'][name] = dict(sha256=sha(p), bytes=p.stat().st_size)
dump(ROOT / 'BENCHMARK_LOCK.json', new)
pred = read(A / 'FROZEN_PREDICTIONS.json')
pred.update(frozen_utc=now(), original_prediction_freeze_utc=previous['frozen_utc'], parent_prediction_freeze_sha256=sha(A / 'FROZEN_PREDICTIONS.json'), benchmark_sha256=sha(ROOT / 'BENCHMARK_LOCK.json'), rebinding_reason='Verifier-only R8 correction after observed release1 failure. All2988 prediction bytes remain unchanged and were frozen before the original measurement; no new blind-test claim.')
dump(ROOT / 'FROZEN_PREDICTIONS.json', pred)
dump(PAYLOAD / 'predictions/FROZEN_PREDICTIONS.json', pred)
write_anchor(sha(ROOT / 'BENCHMARK_LOCK.json'), sha(ROOT / 'FROZEN_PREDICTIONS.json'))
(ROOT / 'RELEASE_ID.txt').write_text('DentalGenCAD-Bench v3 — R8 verifier correction\nbenchmark ' + sha(ROOT / 'BENCHMARK_LOCK.json') + '\npredictions ' + sha(ROOT / 'FROZEN_PREDICTIONS.json') + '\ntrusted_anchor ' + sha(ROOT / 'code/release_anchor.py') + '\nparent_benchmark ' + sha(A / 'BENCHMARK_LOCK.json') + '\n')
dump(ROOT / 'raw/R8_REBINDING.json', dict(utc=now(), changed_operative_files=changed, unchanged_prediction_files=len(previous['files']), unchanged_generator_and_public_inputs=True, parent_benchmark_sha256=sha(A / 'BENCHMARK_LOCK.json'), benchmark_sha256=sha(ROOT / 'BENCHMARK_LOCK.json'), parent_prediction_freeze_sha256=sha(A / 'FROZEN_PREDICTIONS.json'), prediction_freeze_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), limits='Engineering reanalysis of the same observed cohort; not a new preregistered generalization experiment'))
D = PAYLOAD.resolve().parent
for name in changed + ['PREREG_R8.json', 'PREREG_R8.sha256', 'BENCHMARK_LOCK.json', 'FROZEN_PREDICTIONS.json', 'RELEASE_ID.txt', 'code/release_anchor.py']:
    dest = D / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / name, dest)
shutil.copytree(A, D / 'revisions/RELEASE_1', dirs_exist_ok=True)
state('R8_RESEALED_EXPLICITLY', 'Release1 copy failure and numeric-range failure preserved; same prediction bytes rebound to patched verifier', 'Fresh full original-path score and isolated relocated full regeneration/replay')
print('Rebound', len(previous['files']), 'unchanged files; changed', changed)
