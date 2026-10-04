import json, time
from pathlib import Path
import numpy as np
from metrology import write_json, sha
from compare_reference import site_hash, compare
from run_demo import check_frozen
R = Path(__file__).resolve().parents[1]

def run():
    start = time.perf_counter()
    for f in ['PREREG_R3.json', 'DECOMPOSITION_R3.json', 'FROZEN_PREDICTIONS_R3.json']:
        check_frozen(f)
    r1 = json.loads((R / 'R1_RESULTS.json').read_text())
    r2 = json.loads((R / 'R2_RESULTS.json').read_text())
    path = Path(r1['replay_output']) / 'D1/known/occlusal/point_fields.npz'
    oracle = np.load(path)
    truth = oracle['normal_deviation_mm'][0] * 1000
    ref = dict(kind='our_own_fixture', specimen_id='SYNTHETIC_D1', site_sha256=site_hash(oracle), units='um', observable='signed_normal_deviation_at_frozen_design_sites', locator=str(path), bound_locator='PREREG_R3.json: chosen finite-site synthetic error certificate', reference_deviation_um=(truth + 1).tolist(), pointwise_error_bound_um=np.full(len(truth), 2.0).tolist())
    write_json(R / 'raw/REFERENCE_R3.json', ref)
    actual = Path(r2['replay_output']) / 'D1/datum/point_fields.npz'
    report = compare(actual, R / 'raw/REFERENCE_R3.json', 'SYNTHETIC_D1', R / 'raw/REFERENCE_COMPARISON_R3.json')
    lo = np.array(report['manufacture_lower_um'])
    hi = np.array(report['manufacture_upper_um'])
    fails = int(np.sum((truth < lo) | (truth > hi)))
    faulty = []
    for (name, change) in [('unknown_source', {'kind': 'UNKNOWN'}), ('wrong_site', {'site_sha256': 'bad'}), ('wrong_specimen', {'specimen_id': 'OTHER'}), ('missing_bound_locator', {'bound_locator': None}), ('negative_bound', {'pointwise_error_bound_um': [-1.0] * len(truth)})]:
        broken = dict(ref)
        broken.update(change)
        write_json(R / 'raw' / ('bad_reference_' + name + '.json'), broken)
        try:
            compare(actual, R / 'raw' / ('bad_reference_' + name + '.json'), 'SYNTHETIC_D1', R / 'raw/unused_bad_reference.json')
            rejected = False
        except ValueError:
            rejected = True
        faulty.append(dict(name=name, rejected=rejected))
    bad_r = truth + 5
    invalid_enclosure = int(np.sum((truth < bad_r - 2) | (truth > bad_r + 2)))
    gates = dict(true_field_enclosure=fails == 0, interval_width=bool(np.max(abs(hi - lo - 4)) < 1e-10), reject_wrong_sites_specimen_bounds=all((f['rejected'] for f in faulty)), bad_certificate_refuted=invalid_enclosure > 0, synthetic_reference_keeps_unknown=report['empirical_calibration_status'] == 'UNKNOWN_SIMULATED_REFERENCE')
    result = dict(claim_type='capability', gates=gates, finite_site_enclosure_failures=fails, bad_certificate_enclosure_failures=invalid_enclosure, sites=len(truth), faults=faulty, regions=report['regions'], external_referent=dict(kind='our_own_fixture', locator=str(path), compared_quantity='Known simultaneous signed reference field; not independent empirical reference', refutes_us=True), empirical_status='UNKNOWN_NO_PHYSICAL_REFERENCE', prereg_sha256=sha(R / 'PREREG_R3.json'), full_cost={'executed_seconds': time.perf_counter() - start, 'reference_acquisition': 'UNKNOWN'}, precision_control='Direct pointwise interval endpoint arithmetic; exact same information')
    write_json(R / 'R3_RESULTS.json', result)
    combined = json.loads((R / 'results.json').read_text())
    combined['r3'] = result
    write_json(R / 'results.json', combined)
    write_json(R / 'CURRENT_WORK_STATE.json', dict(status='R3_COMPLETE', latest_gate=gates, next_operation='Acquire independent calibrated same-part pointwise reference and external datum bounds; then run compare_reference.py', physical_measurements_acquired=0))
    print('R3 gates', gates, flush=True)
if __name__ == '__main__':
    run()
