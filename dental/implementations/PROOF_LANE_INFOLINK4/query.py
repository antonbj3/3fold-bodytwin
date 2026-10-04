"""Exact source-condition lookup. This returns observations, not fitted predictions."""
from pathlib import Path
import json, sys
H = Path(__file__).resolve().parent

def lookup(registry, link_id, key):
    rows = [r for r in registry.get(link_id, []) if r.get('condition_key') == key]
    if not rows:
        return {'status': 'UNKNOWN_SOURCE_CONDITION', 'link_id': link_id, 'requested_key': key}
    if len(rows) != 1:
        raise ValueError('Ambiguous source key')
    r = rows[0]
    return {k: r[k] for k in ['original_value', 'reported_sd', 'unit', 'n', 'statistic', 'condition_key', 'source_locator', 'source_sha256', 'resolution_level', 'uncertainty_kind']} | {'status': 'SOURCE_OBSERVATION', 'prediction': False}
if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: python3 query.py L01 \'{"material":"NextDent C&B","orientation_deg":0,"postcure_min":120}\'')
    print(json.dumps(lookup(json.loads((H / 'FACIT.json').read_text()), sys.argv[1], json.loads(sys.argv[2])), indent=2, ensure_ascii=False))
