#!/usr/bin/env python3
"""Private graph-use receipts. No job launch, native graph mutation, or evidence fusion."""
import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

WORKSPACES = Path('local_path')
REVIEW_REGISTRY = Path(__file__).resolve().parent / 'GRAPH_USE_AUDIT_20260923/COORDINATOR_REVIEWS.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
        f.write('\n')


def stamp():
    return datetime.now(timezone.utc).isoformat()


def safe_lane(lane):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,100}', lane):
        raise ValueError('Invalid lane ID')
    return lane


def result_path(root, value):
    path = (root / value).resolve()
    if not path.is_relative_to((root / 'results').resolve()) or not path.is_file():
        raise ValueError('Expected existing file under this workspace results/')
    return path


def packet(root, target):
    raw = subprocess.check_output([str(root / 'graph'), 'working', 'packet', '--id', target], text=True)
    return json.loads(raw)


def missing_experiment_fields(p):
    missing = []
    if p.get('numerical_dispatch_allowed') is False or p.get('rank', {}).get('numerical_dispatch_allowed') is False:
        missing.append('numerical_dispatch_not_allowed')
    if p.get('rank', {}).get('blocked_prerequisites'):
        missing.append('unresolved_prerequisites')
    for key in ('best_baseline', 'expected_observation', 'falsifier'):
        if str(p.get(key, '')).strip().upper() in ('', 'NONE', 'NULL', 'UNKNOWN', '—'):
            missing.append(key)
    if not p.get('source_file') and not (p.get('bindings') or p.get('result_bindings')):
        missing.append('source_or_result_binding')
    return missing


def check_feedback(root, record, dispatch=None):
    required = ('target_id', 'result_file', 'sha256', 'review_state', 'outcome',
                'measured_quantity', 'units', 'uncertainty', 'population_regime',
                'preregistered_gate', 'baseline', 'negative_result')
    if any(k not in record for k in required):
        raise ValueError('Incomplete feedback record')
    if dispatch and record['target_id'] != dispatch['target_id']:
        raise ValueError('Feedback target differs from dispatched target')
    source = result_path(root, record['result_file'])
    if digest(source) != record['sha256']:
        raise ValueError('Result hash mismatch')
    reviewed_failure = False
    if record['review_state'] == 'REVIEWED':
        review = result_path(root, record.get('review_artifact', ''))
        if digest(review) != record.get('review_sha256'):
            raise ValueError('Review hash mismatch')
        proof = json.loads(review.read_text())
        if not record.get('review_scope') or any((
            proof.get('target_id') != record['target_id'],
            proof.get('result_sha256') != record['sha256'],
            proof.get('review_scope') != record['review_scope'],
        )):
            raise ValueError('Review scope/result/target mismatch')
        registry = json.loads(REVIEW_REGISTRY.read_text()) if REVIEW_REGISTRY.is_file() else {'reviews': []}
        identity = {'workspace': root.name, 'target_id': record['target_id'],
                    'result_sha256': record['sha256'], 'review_scope': record['review_scope'],
                    'review_sha256': record['review_sha256']}
        if not any(all(row.get(k) == v for k, v in identity.items()) for row in registry['reviews']):
            raise ValueError('Review not admitted by coordinator; retain as pending independent review')
        reviewed_failure = (record.get('reviewer_decision') == 'ACCEPT'
                            and record['outcome'] == 'FAILED'
                            and record['negative_result'] is True
                            and proof.get('decision') == 'ACCEPT_FAILED_GATE')
    elif record['review_state'] != 'PENDING_INDEPENDENT_REVIEW':
        raise ValueError('Unknown review state')
    return {'reviewed_failure_for_task_choice': reviewed_failure,
            'scientific_admission': False,
            'review_authentication': 'coordinator registry plus artifact binding; not cryptographic actor authentication'}


PREVIOUS_INLINE_LIMIT = 50
PREVIOUS_RECENT = 20


def compact_previous(previous):
    """Keep receipts O(1) in target history size.

    Advice is computed from the full history before this call. Reviewed failures
    and integrity issues stay inline; the rest is summarized with an index hash.
    Every FEEDBACK file is preserved, so `graph history --id <target>` recomputes
    the complete list.
    """
    if len(previous) <= PREVIOUS_INLINE_LIMIT:
        return previous, None
    keep = {i for i, x in enumerate(previous)
            if x.get('reviewed_failure_for_task_choice') or x.get('integrity_issue')}
    keep.update(range(len(previous) - PREVIOUS_RECENT, len(previous)))
    outcomes = {}
    for x in previous:
        key = str(x.get('outcome'))[:200]
        outcomes[key] = outcomes.get(key, 0) + 1
    index = hashlib.sha256(''.join(f"{x['receipt']} {x['sha256']}\n" for x in previous).encode()).hexdigest()
    summary = {'schema': 'graph-previous-feedback-summary-v1', 'count': len(previous),
               'inline': len(keep), 'omitted': len(previous) - len(keep),
               'reviewed_failures': sum(bool(x.get('reviewed_failure_for_task_choice')) for x in previous),
               'integrity_issues': sum(bool(x.get('integrity_issue')) for x in previous),
               'outcome_counts': dict(sorted(outcomes.items(), key=lambda kv: -kv[1])[:50]),
               'index_sha256': index,
               'index_rule': 'sha256 over "<receipt> <sha256>\\n" in history order',
               'recompute': 'graph history --id <target_id>'}
    return [previous[i] for i in sorted(keep)], summary


def producer_gain(record):
    """Short verbatim producer gain statement for receipts (unaudited); empty when absent."""
    gain = (record.get('producer_claims') or {}).get('actual_gain')
    return {'producer_gain': gain[:240]} if isinstance(gain, str) and gain else {}


def history(root, target=None):
    rows = []
    for path in sorted((root / 'tasks/graph_runs').glob('*/FEEDBACK_*.json')):
        row = json.loads(path.read_text())
        if target is None or row['record']['target_id'] == target:
            rows.append((path, row))
    return sorted(rows, key=lambda x: x[1]['recorded_at'])


def dispatch(root, target, lane, kind, reason, task_file=None):
    run = root / 'tasks/graph_runs' / safe_lane(lane)
    if run.exists():
        raise ValueError('Lane already has a graph receipt; choose a new run ID')
    p = packet(root, target)
    if kind == 'experiment' and missing_experiment_fields(p):
        raise ValueError('Define missing experiment fields first: ' + ', '.join(missing_experiment_fields(p)))
    previous = []
    for path, row in history(root, target):
        try:
            state = check_feedback(root, row['record'])
        except ValueError as error:
            # A changed historical report is not admissible evidence. Keep the
            # recorded outcome visible and allow only definition/review to repair
            # it; numerical experiment dispatch still fails closed.
            if str(error) != 'Result hash mismatch' or kind == 'experiment':
                raise
            state = {'reviewed_failure_for_task_choice': False,
                     'scientific_admission': False,
                     'integrity_issue': 'RESULT_HASH_MISMATCH',
                     'recorded_outcome_is_unverified': True}
        previous.append({'receipt': str(path), 'sha256': digest(path), **state,
                         'scope': row['record'].get('review_scope'),
                         'outcome': row['record']['outcome'],
                         **producer_gain(row['record'])})
    advice = ('Diagnose the recorded failed gate and choose a discriminating replacement test.'
              if any(x['reviewed_failure_for_task_choice'] for x in previous)
              else p['rank']['action'])
    if any(x.get('integrity_issue') for x in previous):
        advice = 'Historical result integrity is unresolved. Preserve the old recorded outcome as unverified, repair its binding, and do not use it as admitted evidence. ' + advice
    task = None
    if task_file:
        f = task_file.resolve(strict=True)
        task = {'file': str(f), 'sha256': digest(f)}
    previous, summary = compact_previous(previous)
    receipt = {'schema': 'graph-dispatch-v1', 'created_at': stamp(), 'workspace': root.name,
               'lane': lane, 'target_id': target, 'task_kind': kind, 'selection_reason': reason,
               'packet': p, 'previous_feedback': previous, 'next_action_advice': advice,
               'task_file': task, 'launched': False, 'scientific_admission': False,
               'measurement_fields': ['wall_seconds', 'input_tokens_if_observed', 'context_bytes',
                                      'source_citation_errors', 'repeated_refuted_branches',
                                      'baseline_parity', 'falsifiable_next_test', 'scope_errors']}
    if summary:
        receipt['previous_feedback_summary'] = summary
    write_new(run / 'DISPATCH.json', receipt)
    return {'receipt': str(run / 'DISPATCH.json'), 'sha256': digest(run / 'DISPATCH.json'),
            'next_action_advice': advice, 'launched': False}


def feedback(root, lane, record):
    run = root / 'tasks/graph_runs' / safe_lane(lane)
    receipt = json.loads((run / 'DISPATCH.json').read_text())
    checks = check_feedback(root, record, receipt)
    item = {'schema': 'graph-feedback-v1', 'recorded_at': stamp(), 'lane': lane,
            'dispatch_sha256': digest(run / 'DISPATCH.json'), 'record': record, 'checks': checks}
    data = json.dumps(item, sort_keys=True).encode()
    path = run / ('FEEDBACK_' + hashlib.sha256(data).hexdigest()[:16] + '.json')
    write_new(path, item)
    return {'receipt': str(path), **checks}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('project', choices=('dental', 'bodytwin'))
    sub = p.add_subparsers(dest='command', required=True)
    d = sub.add_parser('dispatch')
    d.add_argument('--id', required=True); d.add_argument('--lane', required=True)
    d.add_argument('--kind', choices=('experiment', 'define', 'review'), required=True)
    d.add_argument('--reason', required=True); d.add_argument('--task-file', type=Path)
    f = sub.add_parser('feedback'); f.add_argument('--lane', required=True)
    f.add_argument('--record', type=Path, required=True)
    h = sub.add_parser('history'); h.add_argument('--id')
    a = p.parse_args(); root = WORKSPACES / a.project
    if a.command == 'dispatch':
        out = dispatch(root, a.id, a.lane, a.kind, a.reason, a.task_file)
    elif a.command == 'feedback':
        out = feedback(root, a.lane, json.loads(a.record.read_text()))
    else:
        out = [{'receipt': str(path), **row} for path, row in history(root, a.id)]
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
