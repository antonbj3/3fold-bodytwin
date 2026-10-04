"""Generate proposal locally; install only the explicitly authorized new notes file."""
import argparse
from common import *
NAME = 'constraint_net_dental_r4_20261003.jsonl'

def build():
    n = read('CONSTRAINT_NET_DENTAL_R4.json')
    result_hash = sha(H / 'results.json')
    rows = []
    for (i, e) in enumerate(n['edges']):
        rows.append({'id': 'PROOF_LANE-R4:' + e['id'], 'kind': 'constraint_reference_and_range_proposal', 'source_namespace': 'dental-research/PROOF_LANE_CONSTRAINT_NET_R4', 'edge_id': e['id'], 'between': e['between'], 'inherited_status': e['status'], 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim_type': 'information_link' if 'r4_reference_binding' in e else 'capability', 'depends_on': [], 'binding_semantics': 'proposed context/reference, not admitted physical dependency', 'resolution_level': e['resolution_level'], 'timescale': e['timescale'], 'consumer_decisions': e['decision_ids'], 'r3_contract': {'source': entry('r3_net')['original_path'], 'sha256': entry('r3_net')['sha256'], 'pointer': '/edges/' + str(i)}, 'r4_result': {'file': 'results/PROOF_LANE_CONSTRAINT_NET_R4/results.json', 'sha256': result_hash}, 'range_contract': {'file': 'results/PROOF_LANE_CONSTRAINT_NET_R4/RANGE_CONTRACTS.json', 'sha256': sha(H / 'RANGE_CONTRACTS.json'), 'variables': {k: n['variables'][k]['plausible_range']['status'] for k in e['between']}}, 'new_primary_reference': e.get('r4_reference_binding'), 'physical_status_changed': False, 'rigorous_enclosure': 'No new physical/model transfer enclosure; retain R3 limitations.'})
    (H / NAME).write_text(''.join((json.dumps(x, ensure_ascii=False) + '\n' for x in rows)))
    return rows

def install():
    dst = H.parent.parent / 'notes/expansion' / NAME
    data = (H / NAME).read_bytes()
    if dst.exists() and dst.read_bytes() != data:
        raise RuntimeError('NEW_BINDING_ALREADY_EXISTS_WITH_DIFFERENT_BYTES: preserve and resolve explicitly')
    if not dst.exists():
        with dst.open('xb') as f:
            f.write(data)
    write('BINDINGS_RECEIPT.json', {'file': str(dst), 'sha256': sha(dst), 'rows': len((H / NAME).read_text().splitlines()), 'new_file_only': True, 'source_graphs_modified': False})
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--install', action='store_true')
    a = ap.parse_args()
    build()
    if a.install:
        install()
