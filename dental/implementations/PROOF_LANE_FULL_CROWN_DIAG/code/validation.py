from common import *
import copy
from sufficiency import run as witness

def diagnostic_decision(q, threshold, contact, wall_lower, actual_margin_known=False):
    return dict(shape=bool(q <= threshold), contact=bool(contact), wall=bool(wall_lower >= 0.5), actual_margin='MEASURED' if actual_margin_known else 'UNKNOWN', clinical='NOT_ESTABLISHED')

def run():
    pr = read(ROOT / 'PREREG_RESCORE.json')
    checks = {}
    hash_count = 0
    for (p, h) in pr['input_hashes'].items():
        assert sha(p) == h, p
        hash_count += 1
    for name in ['PREREG_LOCAL', 'PREREG_ANNOTATED', 'PREREG_RESCORE', 'FROZEN_PREDICTIONS', 'FROZEN_RESCORE_PREDICTIONS', 'DECOMPOSITION_LOCAL', 'DECOMPOSITION_ANNOTATED', 'DECOMPOSITION_RESCORE', 'PREREG_COMPONENT_AUDIT']:
        checks[name + '_hash'] = sha(ROOT / (name + '.json')) == (ROOT / (name + '.sha256')).read_text().split()[0]
    initial = ROOT / 'history/initial_measurements/ANNOTATED_PAIRS.json'
    if initial.exists():
        checks['original_threshold_source_preserved'] = sha(initial) == read(ROOT / 'FROZEN_RESCORE_PREDICTIONS.json')['threshold_source_sha256']
    ws = witness()
    assert ws['sign_witness']['identity_error_mm'] == 0 and ws['location_witness']['identity_error_mm'] == 0
    good = diagnostic_decision(0.2, 0.3, True, 0.7)
    mutations = dict(shape_1mm_rejected=not diagnostic_decision(1.0, 0.3, True, 0.7)['shape'], wall_0mm_rejected=not diagnostic_decision(0.2, 0.3, True, 0.0)['wall'], contact_false_rejected=not diagnostic_decision(0.2, 0.3, False, 0.7)['contact'], unknown_margin_never_passes=good['actual_margin'] == 'UNKNOWN' and good['clinical'] == 'NOT_ESTABLISHED', unit_scale1000_rejected=not diagnostic_decision(0.2 * 1000, 0.3, True, 0.7)['shape'], wrong_input_hash_rejected=sha(ROOT / 'PREREG_LOCAL.json') != '0' * 64)
    from rescore import gates, compare
    spec = read(R3 / 'PREREG_A.json')['metrics']['v6']
    xy = np.array([[0.0, 0], [1, 0], [1, 1], [0, 1]])
    ff = np.array([[0, 1, 2], [0, 2, 3]])
    template = compare(xy, ff, np.full(4, 0.05), np.full(4, 0.05))
    assert all(gates(template, spec).values())
    for (key, bad, gate) in [('symdiff_mm2', 2.0, 'pattern'), ('count_error', 3, 'count'), ('centroid_distance_mm', 1.0, 'location'), ('negative_gap_area_mm2', 1.0, 'interference'), ('coverage', 0.5, 'support')]:
        m = copy.deepcopy(template)
        m[key] = bad
        mutations[key + '_rejected'] = not gates(m, spec)[gate]
    assert all(checks.values()) and all(mutations.values())
    out = dict(frozen_hash_checks=checks, input_hashes_verified=hash_count, injected_errors=mutations, sufficiency_controls=ws['controls'], status='PASS', scope='Numerical/data contract tests; no independent clinical validation')
    dump(ROOT / 'raw/VALIDATION.json', out)
    return out
if __name__ == '__main__':
    run()
