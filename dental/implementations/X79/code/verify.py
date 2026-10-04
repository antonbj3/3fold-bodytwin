from common import *
from experiment import metrics
from report_parser import parse
from data import load_geometry
import zipfile, numpy as np

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    checks = {}
    for k in ['R1', 'R3']:
        frozen = read('FROZEN_PREDICTIONS_' + k + '.json')
        assert sha(frozen['predictions_path']) == frozen['predictions_sha256']
        assert sha(frozen['model_path']) == frozen['model_sha256']
        rows = read('raw/' + k + '_TEST_ROWS.json')
        bad = metrics(rows, [1 - r['y'] for r in rows])
        checks[k + '_wrong_value_rejected'] = bad['balanced_accuracy'] == 0 and bad['sensitivity'] == 0 and (bad['specificity'] == 0)
        checks[k + '_frozen_prediction_hash'] = True
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    sets = [set(m[k]) for k in ['train', 'calibration', 'test']]
    checks['patient_disjoint'] = not any((sets[i] & sets[j] for i in range(3) for j in range(i)))
    with zipfile.ZipFile(ZIP) as z:
        checks['source_hashes'] = True
        for key in ['R1', 'R3']:
            for r in read('raw/' + key + '_TEST_ROWS.json'):
                assert digest(z.read(r['source_member'])) == r['source_sha256']
                assert r['source_member'].split('/')[0] == r['case_id']
    checks['material_port_wrong_identity'] = all(read('raw/R2_RESULTS.json')['injections'].values())
    checks['same_information_material_control'] = read('raw/R2_RESULTS.json')['matched_control_parity_mismatches'] == 0
    checks['same_information_material_wrong_value_rejected'] = read('raw/R2_RESULTS.json')['injections']['same_information_control_injected_state_rejected']
    checks['R5_native_and_site_scope_injections'] = all(read('raw/R5_RESULTS.json')['injections'].values())
    checks['silent_report_unknown'] = parse('The dental arches are described.')['assertions'] == []
    checks['negation_not_inherited_across_but'] = parse('No restorations are present, but there is crowding.')['assertions'][0]['polarity'] == 0
    checks['primary_not_permanent'] = parse('A restoration is present on tooth75.')['assertions'] == []
    checks['hedge_refused'] = parse('A restoration may be present on tooth16.')['assertions'] == []
    checks['sealant_not_restoration'] = parse('Pit-and-fissure sealants are present on16.')['assertions'] == []
    checks['sealant_or_restoration_refused'] = parse('Sealants/conservative restorations are present on16.')['assertions'] == []
    checks['global_other_negative_refused'] = parse('No other restorative treatments are present.')['assertions'] == []
    checks['source_gold_precision'] = read('raw/PARSER_GOLD_VALIDATION.json')['precision'] >= 0.95
    checks['sufficiency_exact_identity'] = read('raw/SUMMARY_SUFFICIENCY.json')['pass_gate']
    if (ROOT / 'raw/R4_RESULTS.json').exists():
        checks['pointset_distance_control'] = read('raw/R4_RESULTS.json')['pointset_comparator']['pass_gate']
    write('raw/VERIFICATION.json', dict(pass_gate=all(checks.values()), checks=checks, cost=cost(st, cpu), scope='Computation/source alignment; not independent clinical or physical validation'))
    assert all(checks.values()), checks
    print('Verification checks', len(checks), 'PASS', flush=True)
if __name__ == '__main__':
    main()
