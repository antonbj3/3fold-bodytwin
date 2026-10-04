"""Delivery integrity only; no new scientific threshold or independent review."""
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parsed = {p.name: json.loads(p.read_text()) for p in ROOT.glob('*.json')}
    frozen = {}
    for p in list(ROOT.glob('PREREG*.sha256')) + list(ROOT.glob('FROZEN_PREDICTIONS*.sha256')):
        data = p.with_suffix('.json')
        frozen[data.name] = sha(data) == p.read_text().strip()
    assert frozen and all(frozen.values()), 'Frozen preregistration/prediction drift'
    source_checks = {}
    for name in ['SOURCE_MANIFEST.json', 'SOURCE_MANIFEST_R3.json']:
        for row in parsed[name]:
            source_checks[row['local']] = sha(ROOT / row['local']) == row['sha256']
    assert all(source_checks.values()), 'Local input drift'
    with (ROOT / 'DECIDABILITY_TABLE.csv').open() as f:
        rows = list(csv.DictReader(f))
    parents = [r for r in rows if r['decision'].split(':', 1)[-1].startswith(('physical_', 'absolute_', 'biological_', 'cement_', 'bone_'))]
    canal = [r for r in rows if r['decision'].startswith('X5:') and r not in parents]
    reference = [r for r in rows if r['decision'].startswith('REFERENCE:')]
    results = parsed['results.json']
    counts = Counter((r['class'] for r in canal))
    assert len(parents) == results['parent_decisions']['covered'] == 8
    assert all((r['class'] == '2' and r['binding_quantity'] for r in parents))
    assert len(canal) == results['canal_conditional']['query_instances'] == 4132
    assert counts['1'] == results['canal_conditional']['class1'] == 3648
    assert counts['2'] == results['canal_conditional']['class2'] == 484
    assert len(reference) == 1 and reference[0]['class'] == '3'
    assert len(rows) == 4141 and all((r['physical_status'] == 'UNKNOWN' for r in rows))
    assert all((r['conditional_floor_resolution_level'] == 'PHENOMENOLOGICAL' and r['numerical_pitch_resolution_level'] == 'PER_POINT' and (r['signed_margin_resolution_level'] == 'PER_SURFACE_REGION') for r in canal))
    round_checks = {}
    for name in ['RESULTS_R2.json', 'RESULTS_R3.json', 'RESULTS_R4.json']:
        d = parsed[name]
        round_checks[name] = all(d['gates'].values()) and all(d['mutation_rejections'].values())
    assert all(round_checks.values())
    feedback = parsed['GRAPH_FEEDBACK.json']
    assert feedback['sha256'] == sha(ROOT / 'results.json')
    assert feedback['review_state'] == results['review_state'] == 'PENDING_INDEPENDENT_REVIEW'
    expansion = ROOT.parents[1] / 'notes/expansion/decidability_20261003.jsonl'
    proposed = [json.loads(x) for x in expansion.read_text().splitlines() if x.strip()]
    assert len(proposed) == 8 and len({r['id'] for r in proposed}) == 8
    assert expansion.read_bytes() == (ROOT / 'GRAPH_EXPANSION_PROPOSAL_COPY.jsonl').read_bytes()
    report = {'scope': 'Delivery integrity; neither a new empirical gate nor independent scientific review', 'checked_utc': datetime.now(timezone.utc).isoformat(), 'root_json_files_parsed': len(parsed), 'frozen_hashes': frozen, 'source_hashes': source_checks, 'table_rows': len(rows), 'canal_questions': len(canal), 'conditional_classes': dict(counts), 'parent_class2': len(parents), 'analytical_reference_class3': len(reference), 'physical_status': 'UNKNOWN', 'round_gates_and_injections': round_checks, 'graph_binding_matches': True, 'review_state': feedback['review_state'], 'expansion_proposals': len(proposed), 'expansion_sha256': sha(expansion), 'scientific_admission': False}
    (ROOT / 'FINAL_DELIVERY_AUDIT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'table_rows': len(rows), 'conditional_classes': dict(counts), 'physical_status': 'UNKNOWN', 'graph_binding_matches': True}))
if __name__ == '__main__':
    main()
