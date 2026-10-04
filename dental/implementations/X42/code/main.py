"""One-command locked replay, researcher report and figure."""
import argparse, json, sys
from common import *
from freeze_predictions import verify_generator
from referee_v2 import eval_once

def run():
    ap = argparse.ArgumentParser()
    ap.add_argument('--first-measurement', action='store_true')
    ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--dev-demo-only', action='store_true')
    a = ap.parse_args()
    verify_generator()
    wrappers = json.load(open(ROOT / 'FROZEN_EXECUTION_WRAPPERS.json'))
    if sha(ROOT / 'FROZEN_EXECUTION_WRAPPERS.json') != (ROOT / 'FROZEN_EXECUTION_WRAPPERS.sha256').read_text().split()[0]:
        raise ValueError('execution wrapper anchor mismatch')
    for (name, entry) in wrappers['files'].items():
        if sha(ROOT / name) != entry['sha256']:
            raise ValueError('execution wrapper changed: ' + name)
    pred = json.load(open(ROOT / 'FROZEN_PREDICTIONS.json'))
    receipt = pred['manifest_sha256']
    if sha(Path(pred['manifest_path'])) != receipt:
        raise ValueError('prediction manifest changed')
    if a.verify_only:
        from integrity import Bundle, predictions
        from release_anchor import BENCHMARK_SHA256
        Bundle(BENCH)
        predictions(DATA / 'PREDICTIONS', receipt, [r['key'] for r in split()['metadata']], BENCHMARK_SHA256)
        print('locked generator, benchmark and predictions verified')
        return
    if a.dev_demo_only:
        from demo_generate import run_demo
        run_demo()
        return
    replay = not a.first_measurement
    if not replay and (ROOT / 'raw/QUERY_LEDGER.json').exists():
        raise ValueError('first query already charged; use locked replay')
    agg = eval_once(DATA / 'PREDICTIONS', receipt, replay=replay)
    if replay:
        first = json.load(open(ROOT / 'raw/FINAL_AGGREGATES.json'))
        from util import digest

        def stable(o):
            if isinstance(o, dict):
                return {k: stable(v) for (k, v) in o.items() if k not in ['seconds', 'peak_rss_MiB', 'replay_of_same_frozen_predictions']}
            if isinstance(o, list):
                return [stable(v) for v in o]
            return o
        if digest(stable(agg)) != digest(stable(first)):
            raise ValueError('replay changed scientific aggregates')
        dump(ROOT / 'raw/REPLAY_RECEIPT.json', dict(PASS=True, scientific_digest=digest(stable(agg)), prediction_freeze_sha256=receipt, new_hidden_queries=0))
        agg = first
    from demo_generate import run_demo
    run_demo()
    from x42_report import report
    report(agg)
    if replay:
        print('Same frozen prediction replay; no new hidden query.')
    print('X42 complete digital evaluation; results.json, README_DEMO.md, figures/x42.png. PENDING_INDEPENDENT_REVIEW.')
if __name__ == '__main__':
    run()
