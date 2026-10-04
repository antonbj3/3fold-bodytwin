import sys, json
from common import *

def main():
    demo_only = '--demo-only' in sys.argv
    rs = read('results.json')
    cs = read('SOURCE_CONTRACTS.json')
    regs = read('OBSERVATION_REGISTRY.json')
    r3 = read('R3_RESULTS.json')
    fault = read('FAULT_INJECTION.json')
    for p in ROOT.glob('PREREG*.json'):
        assert sha(p) == (ROOT / (p.name + '.sha256')).read_text().strip()
    for p in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        assert sha(p) == (ROOT / (p.name + '.sha256')).read_text().strip()
    for (path, rec) in read('SOURCE_MANIFEST.json').items():
        assert sha(path) == rec['sha256']
    lock = ROOT / 'REPRODUCTION_INPUTS_LOCK.json'
    assert sha(lock) == (ROOT / (lock.name + '.sha256')).read_text().strip()
    for (path, rec) in read(lock.name)['inputs'].items():
        assert sha(path) == rec['sha256']
    assert rs['claim_type'] == 'information_link' and rs['established_source_links'] == 4 and (rs['blocked_candidate_links'] == 1)
    assert len(regs) == 51 and fault['all_pass'] and (rs['faults_rejected'] == 377)
    assert r3['independent_point_inside'] and r3['bad_plus1C_rejected']
    assert rs['summary_identity_error_max'] == 0 and rs['full_native_chains_empirically_closed'] == 0
    assert all((v['summary_bitwise_equal'] for v in rs['sufficiency_tests'].values()))
    assert not rs['R1']['L01']['source_gate'] and (not rs['R1']['L04']['coefficient_consistency_gate_pass'])
    assert all((not x['frozen_gate_pass'] for x in rs['R1']['L02']['paired_ratios']))
    assert all((x['bad_value_rejected'] for x in r3['control_fault_tests']))
    assert read('FACIT.json')['all_pass']
    for p in ['RESULTS.md', 'README_DEMO.md', 'figure.png', 'figure.svg', 'figure.pdf', 'information_links.preview.jsonl']:
        assert (ROOT / p).stat().st_size > 0, p
    nodes = [json.loads(l) for l in (ROOT / 'information_links.preview.jsonl').read_text().splitlines()]
    assert len(nodes) == 5
    assert all((n['claim_type'] == 'information_link' and n['review_state'] == 'PENDING_INDEPENDENT_REVIEW' for n in nodes))
    assert sum((n['source_measurement_gate'] for n in nodes)) == 4
    if not demo_only:
        install = read('GRAPH_INSTALL_CHECK.json')
        assert install['all_imported'] and (not install['native_graph_files_changed'])
        assert sha(install['installed_path']) == install['sha256']
        assert Path(install['installed_path']).read_bytes() == (ROOT / 'information_links.preview.jsonl').read_bytes()
    bytes_total = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    assert bytes_total < 3000000000
    put('DELIVERY_CHECK.json', {'all_pass': True, 'demo_only': demo_only, 'quantity_records': 51, 'primary_links': 4, 'blocked_candidates': 1, 'corruptions_rejected': 377, 'summary_identity_error_max': 0, 'artifact_bytes': bytes_total, 'large_arrays': [], 'result_sha256': sha(ROOT / 'results.json'), 'graph_install_checked': not demo_only, 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print('Delivery validation PASS;377 refutable controls,4 primary links,1 blocked; no full-chain closure.')
if __name__ == '__main__':
    main()
