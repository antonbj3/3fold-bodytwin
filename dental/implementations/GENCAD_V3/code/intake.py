"""Trusted referee intake. One external test query total per release.

Participant code never receives this directory or its receipts. Its generator
identity is recorded; changing it does not reset the release query budget.
"""
import re, os, fcntl, shutil, json
from pathlib import Path
from integrity import Bundle, IntegrityError, safe_bytes, h, file_set, load_npz
from scorer import decode, score_submission
from util import dump, now, sha

def submit(root, candidate, generator_sha256):
    if not re.fullmatch('[0-9a-f]{64}', generator_sha256):
        raise IntegrityError('frozen generator SHA256 required')
    b = Bundle(root)
    from release_anchor import BENCHMARK_SHA256
    cases = {Path(p).stem: ts for (p, ts) in b.tasks()}
    candidate = Path(candidate)
    if file_set(candidate) != {key + '.npz' for key in cases}:
        raise IntegrityError('submission must cover exactly every case')
    base = b.payload / 'referee_intake'
    base.mkdir(mode=448, exist_ok=True)
    with (base / 'ledger.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        ledger = base / 'QUERY_LEDGER.json'
        old = json.loads(ledger.read_text()) if ledger.exists() else []
        if old:
            raise IntegrityError('one-query external heldout budget exhausted for this release')
        dest = base / generator_sha256
        if dest.exists():
            raise IntegrityError('pending/failed intake exists; referee must inspect, no automatic reset')
        (dest / 'candidate').mkdir(parents=True)
        files = {}
        for (key, ts) in cases.items():
            blob = safe_bytes(candidate, key + '.npz', 2000000)
            decode(load_npz(blob, max_expanded=2000000), ts)
            name = 'candidate/' + key + '.npz'
            (dest / name).write_bytes(blob)
            files[name] = dict(sha256=h(blob), bytes=len(blob))
        manifest = dict(frozen_utc=now(), benchmark_sha256=BENCHMARK_SHA256, generator_sha256=generator_sha256, participants=['candidate'], files=files)
        dump(dest / 'FROZEN_PREDICTIONS.json', manifest)
        receipt = sha(dest / 'FROZEN_PREDICTIONS.json')
        old.append(dict(generator_sha256=generator_sha256, prediction_freeze_sha256=receipt, accepted_utc=now(), test_queries_charged=1))
        dump(ledger, old)
    scores = score_submission(b.root, dest, receipt, base / (generator_sha256 + '_scores'))
    return dict(benchmark_sha256=BENCHMARK_SHA256, receipt_sha256=receipt, summaries=scores['summaries'])
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('candidate_directory')
    ap.add_argument('generator_sha256')
    a = ap.parse_args()
    print(json.dumps(submit(Path(__file__).resolve().parents[1], a.candidate_directory, a.generator_sha256), indent=2))
