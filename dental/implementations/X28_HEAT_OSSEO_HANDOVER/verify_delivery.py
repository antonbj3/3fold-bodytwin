import json
from pathlib import Path
from common import ROOT, load, sha, verify_freeze, dump, state

def main():
    required = ['CURRENT_WORK_STATE.json', 'RESULTS.md', 'results.json', 'COMMANDS.md', 'HANDOFF.md', 'README_DEMO.md', 'GRAPH_FEEDBACK.json', 'MINIMUM_MEASUREMENT_CONTRACT.json', 'SOURCE_MANIFEST.json', 'PROPOSED_INFORMATION_EDGES.json', 'SOURCE_TIMING_CORRECTION.json', 'figures/X28_measured_handover.png', 'figures/X28_K3_thermal_port.png']
    missing = [p for p in required if not (ROOT / p).is_file()]
    if missing:
        raise ValueError('Missing delivery: ' + str(missing))
    frozen = ['FROZEN_PREDICTIONS.json'] + [f'FROZEN_PREDICTIONS_R{i}.json' for i in [2, 3, 4]] + [f'PREREG_R{i}.json' for i in [1, 2, 3, 4]]
    for p in frozen:
        verify_freeze(p)
    (r1, r2, r3, r4) = [load(f'results_R{i}.json') for i in range(1, 5)]
    assert all(r1['faults'].values()) and all(r2['checks'].values()) and all(r3['fault_controls'].values()) and all(r4['fault_controls'].values())
    assert r3['frozen_prediction_check']['all_passed'] is False
    assert r3['frozen_prediction_check']['checks'][0]['expected'] == 39.1
    assert abs(r3['measured_contrasts']['viability_500minus1000_percentage_points'] - 37.59) < 1e-08
    assert all((p['implant_ISQ'] == 'UNKNOWN' and p['healing_delay_weeks'] == 'UNKNOWN' for p in load('BIOLOGICAL_HANDOVER_PORTS.json')))
    assert all((p['predicted_ISQ'] == 'UNKNOWN' for p in load('THERMAL_HANDOVER_STATES.json')))
    feedback = load('GRAPH_FEEDBACK.json')
    assert feedback['sha256'] == sha(ROOT / 'results.json') and feedback['review_state'] == 'PENDING_INDEPENDENT_REVIEW'
    for source in load('SOURCE_MANIFEST.json'):
        assert sha(source['path']) == source['sha256'], source['path']
    total = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    assert total < 3000000000
    manifest = [dict(path=str(p.relative_to(ROOT)), sha256=sha(p), bytes=p.stat().st_size) for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in str(p) and (p.name not in ['DELIVERY_MANIFEST.json', 'DELIVERY_CHECKS.json', 'CURRENT_WORK_STATE.json']) and (not p.name.endswith('.log'))]
    dump('DELIVERY_MANIFEST.json', manifest)
    dump('DELIVERY_CHECKS.json', dict(technical_checks_passed=True, frozen_files_verified=len(frozen), source_hashes_verified=True, feedback_hash_verified=True, all_fault_controls_passed=True, R3_wrong_frozen_contrast_preserved=True, R3_prediction_passed=False, physiological_prediction_validated=False, clinical_safe_protocol_known=False, total_output_bytes=total, scientific_status='PENDING_INDEPENDENT_REVIEW'))
    state('LOCAL_DEMO_VERIFIED_EMPIRICAL_CHAIN_OPEN', 'technical/source/fault checks pass;R3frozencontrastFAIL preserved;healing/ISQUNKNOWN', 'Independent review and matched retained-wall thermometry/finalhistology/serialstability;graph owner refresh then bind')
    print('Delivery verified; R3 frozen contrast remains FAIL; physiological transfer remains UNKNOWN.')
if __name__ == '__main__':
    main()
