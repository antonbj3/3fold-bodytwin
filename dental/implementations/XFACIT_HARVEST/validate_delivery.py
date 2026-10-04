from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib, json, collections, datetime
P = Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
r = json.loads((P / 'results.json').read_text())
obs = list(map(json.loads, (P / 'VERIFIED_OBSERVATIONS.jsonl').read_text().splitlines()))
ver = list(map(json.loads, (P / 'VERIFICATIONS.jsonl').read_text().splitlines()))
draft = list(map(json.loads, (P / 'facit_harvest_20261004.jsonl').read_text().splitlines()))
feedback = json.loads((P / 'GRAPH_FEEDBACK.json').read_text())
receipt = json.loads((P / 'GRAPH_FEEDBACK_RECEIPT.json').read_text())
receipt_record = json.loads(Path(receipt['receipt']).read_text())['record']
assert r['harvest']['selected_reports'] == len(ver) == 52
assert sum((v['source_numerical_adoption'] for v in ver)) == 40
assert len(obs) == len({o['id'] for o in obs}) == 269
assert len(draft) == len({d['id'] for d in draft}) == 606
assert len({d['existing_edge_id'] for d in draft}) == 24
assert all((o['population'] and o['protocol'] and o['unit'] and o['reference']['locator'] for o in obs))
assert all((sha(o['reference']['source_file']) == o['reference']['source_sha256'] for o in obs))
assert all((d['no_scientific_admission'] and (not d['plausible_range']['target_numeric_admission']) for d in draft))
assert all((d['resolution_level'] == d['observation']['resolution_level'] for d in draft))
assert all((d['timescale'] in ['SIMULTANEOUS', 'HANDOVER'] for d in draft))
assert receipt_record == feedback
assert sha(P / 'results.json') == feedback['sha256']
assert sha(P / 'facit_harvest_20261004.jsonl') == feedback['proposal_sha256']
assert sha(_release_expand('@DENTAL_INPUT_ROOT@/workspace/') + feedback['proposal_file']) == feedback['proposal_sha256']
assert json.loads((P / 'REPLAY_VERIFICATION.json').read_text())['passed']
lock = json.loads((P / 'R4_LOCK.json').read_text())
assert sha(lock['source']) == lock['sha256']
static = []
for side in P.glob('*.sha256'):
    if side.name == 'CURRENT_WORK_STATE.json.sha256':
        continue
    expected = side.read_text().split()[0]
    target = P / side.name[:-7]
    assert target.exists() and sha(target) == expected, side.name
    static.append(side.name)
initial = json.loads((P / 'SOURCE_FETCH_LOG.json').read_text())['requests']
additional = json.loads((P / 'SOURCE_FETCH_R3_LOG.json').read_text())
cost = {'known_logged_source_requests': len(initial) + len(additional), 'initial_status': dict(collections.Counter((x['status'] for x in initial))), 'additional_status': dict(collections.Counter((x['status'] for x in additional))), 'request_wall_s_sum': sum((x.get('wall_s', 0) for x in initial + additional)), 'unlogged_inline_source_requests': 'UNKNOWN count/time; primary response snapshots hashed', 'web_search_count': 'UNKNOWN; source acquisition manifests do not account for every research search', 'pre_prereg_instruction_reads_and_search_wall_s': 'UNKNOWN', 'discovery_and_validation_elapsed_from_R1_to_delivery_s': (datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(json.loads((P / 'PREREG_XFACIT_R1.json').read_text())['frozen_at'])).total_seconds(), 'numerical_fit': 0, 'human_questions': 0, 'GPU': False, 'demo_threads': 1, 'demo_runtime': json.loads((P / 'REPLAY_VERIFICATION.json').read_text())['runtime'], 'fallback': 'Primary XML failures, corrected identifiers, publisher HTML, visual PDF;12report candidates retained as no numerical adoption', 'resolution_level': 'PHENOMENOLOGICAL'}
(P / 'COST_LEDGER.json').write_text(json.dumps(cost, indent=2) + '\n')
cache = Path(_release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest'))
files = [p for p in cache.rglob('*') if p.is_file()]
total = sum((p.stat().st_size for p in files))
assert total < 3000000000
(P / 'SOURCE_CACHE_MANIFEST.json').write_text(json.dumps({'base': str(cache), 'total_bytes': total, 'files': [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in sorted(files)]}, indent=2) + '\n')
out = {'passed': True, 'validated_at': now, 'source_values': len(obs), 'adjudications': len(ver), 'draft_facets': len(draft), 'unique_edges': 24, 'source_cache_bytes': total, 'static_frozen_sidecars_checked': static, 'canonical_result_hash': sha(P / 'results.json'), 'graph_feedback_bound': True, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'R4_source_unchanged': True, 'draft_export_identical': True, 'zero_physical_target_admissions': True, 'scope': 'Artifact, provenance and executed-gate audit; independent scientific review remains pending'}
(P / 'FINAL_VALIDATION.json').write_text(json.dumps(out, indent=2) + '\n')
state = {'updated_at': now, 'phase': 'DELIVERED_PENDING_INDEPENDENT_REVIEW', 'latest_gate': r['outcome'], 'claim_type': 'information_link', 'in_progress': None, 'next_operation': 'Coordinator reviews exact source/quantity/protocol bindings; then XFACIT_PROTOCOL_BRIDGE at K20, no automatic graph admission', 'last_validation': 'FINAL_VALIDATION.json', 'result_sha256': sha(P / 'results.json'), 'counts': r['harvest'], 'mapping': r['mapping'], 'open_obstacle': 'No matched physical target protocol/range; individual field uncertainty unresolved', 'feedback_receipt': receipt['receipt'], 'handoff': 'HANDOFF.md'}
(P / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, ensure_ascii=False, indent=2) + '\n')
(P / 'CURRENT_WORK_STATE.json.sha256').write_text(sha(P / 'CURRENT_WORK_STATE.json') + '\n')
print(json.dumps(out))
