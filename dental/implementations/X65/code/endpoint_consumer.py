"""Exact protocol identity for a published benchmark endpoint; no risk inference."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def endpoint_query(query):
    required = ['system_id', 'diameter_mm', 'material', 'connection', 'angle_deg', 'R', 'horizon_cycles', 'quantity', 'unit']
    if any((k not in query for k in required)):
        return {'status': 'UNKNOWN', 'reason': 'missing full configuration/protocol identity'}
    data = json.loads((ROOT / 'raw/published_endpoints.json').read_text())
    row = next((r for r in data['records'] if r['id'] == query['system_id']), None)
    if not row:
        return {'status': 'UNKNOWN', 'reason': 'system has no primary measurement'}
    if not row['strict_eligible']:
        return {'status': 'UNKNOWN', 'reason': 'strict primary protocol not complete', 'rejections': row['strict_reasons']}
    for k in ['diameter_mm', 'material', 'connection', 'angle_deg', 'R', 'horizon_cycles']:
        if query[k] != row[k]:
            return {'status': 'UNKNOWN', 'reason': 'configuration/protocol mismatch: ' + k}
    if query['quantity'] != 'maximum_cyclic_force' or query['unit'] != 'N':
        return {'status': 'UNKNOWN', 'reason': 'endpoint is not risk, life, stress amplitude or another unit'}
    if 'reference_force_N' in query and query['reference_force_N'] != row['force_N']:
        return {'status': 'UNKNOWN', 'reason': 'asserted endpoint contradicts primary cell'}
    c = query.get('force_class_N')
    return {'status': 'PUBLISHED_ENDPOINT', 'maximum_cyclic_force_N': row['force_N'], 'horizon_cycles': row['horizon_cycles'], 'support_at_requested_class': 'AT_OR_BELOW_REPORTED_ENDPOINT' if c is not None and c <= row['force_N'] else 'ABOVE_REPORTED_ENDPOINT_NOT_DEMONSTRATED', 'class_N': c, 'criterion': row['criterion'], 'locator': row['locator'], 'table_locator': row['table_locator'], 'resolution': 'PER_TOOTH', 'risk_certificate': False, 'new_configuration_prediction': False}
