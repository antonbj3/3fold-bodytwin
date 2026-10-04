from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '4'
from pathlib import Path
import json, hashlib, datetime, time
import numpy as np
from model import response, window
R = Path(__file__).resolve().parents[1]

def write(p, x):
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def predict():
    prereg = json.load(open(R / 'PREREG_R1.json'))
    geom = json.load(open(R / 'inputs/ARCH_GEOMETRY.json'))
    intervals = json.load(open(_release_expand('@DENTAL_INPUT_ROOT@/workspace/notes/aligner_intervals/ALIGNER_INTERVALS.json')))
    rows = []
    refs = []
    controls = []
    for a in geom['arches']:
        for k0 in sorted(a['teeth'], key=int):
            k = int(k0)
            p = k % 10
            record = intervals[1 if p <= 3 else 2]
            central = record['value_central']
            bounds = [record['value_min'], record['value_max']]
            for delta in [0.2, 0.25]:
                nom = response(a, k, delta, 2746, 0.75)
                b = response(a, k, delta, 2746, central)
                full = response(a, k, delta, 2189, central)
                if full['status'].startswith('REFUSED'):
                    rows.append({'case': a['case'], 'jaw': a['jaw'], 'active_fdi': k, 'activation_mm': delta, 'status': full['status']})
                    continue
                lo = response(a, k, delta, 2189, bounds[0])
                hi = response(a, k, delta, 2189, bounds[1])
                con = response(a, k, delta, 2189, central, True)
                controls.append(max((np.max(np.abs(np.array(full['wrenches'][t]) - np.array(con['wrenches'][t]))) for t in full['wrenches'])))
                reg = []
                for idx in [3, 4, 5, 6, 7]:
                    rr = intervals[idx]
                    g = response(a, k, delta, 2189, rr['value_central'])
                    reg.append({'region': rr['region'], 'source_doi': rr['source_doi'], 'thickness_mm': rr['value_central'], 'force_N': g['active_force_norm_N'], 'moment_Nmm': g['active_moment_norm_Nmm'], 'statistic': rr['statistic'], 'tooth_by_region_joint_status': 'NOT_AVAILABLE; marginal regional sensitivity only'})
                r = {'case': a['case'], 'dataset': a['dataset'], 'label_status': a['label_status'], 'jaw': a['jaw'], 'tooth_type': a['teeth'][k0]['tooth_type'], 'active_fdi': k, 'activation_mm': delta, 'status': full['status'], 'nominal': nom, 'thickness_only': b, 'formed': full, 'formed_IQR_force_interval_N': [lo['active_force_norm_N'], hi['active_force_norm_N']], 'IQR_is_not_confidence_or_guaranteed_range': True, 'descriptive_bodily_range_screen': window([lo['active_force_norm_N'], hi['active_force_norm_N']]), 'force_change_thickness_only_percent': 100 * (b['active_force_norm_N'] / nom['active_force_norm_N'] - 1), 'force_change_thickness_and_E_percent': 100 * (full['active_force_norm_N'] / nom['active_force_norm_N'] - 1), 'thickness_record_index': 1 if p <= 3 else 2, 'regions': reg, 'patient_force_status': 'UNKNOWN_PHYSICAL_MODEL_DISCREPANCY'}
                rows.append(r)
            refs.append({'case': a['case'], 'fdi': k, 'thickness_record_index': 1 if p <= 3 else 2, 'doi': record['source_doi'], 'locator': record['locator'], 'transferred_not_specimen_matched': True})
    out = {'schema': 'x19-predictions-r1-v1', 'claim_type': 'information_link', 'rows': rows, 'thickness_bindings': refs, 'hermite_max_wrench_difference': max(controls), 'sha256_inputs': {'prereg': sha(R / 'PREREG_R1.json'), 'geometry': sha(R / 'inputs/ARCH_GEOMETRY.json'), 'intervals': sha(Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/notes/aligner_intervals/ALIGNER_INTERVALS.json')))}, 'code_sha256': {p.name: sha(p) for p in [R / 'code/model.py', R / 'code/predict_r1.py']}}
    return out
if __name__ == '__main__':
    t = time.perf_counter()
    p = R / 'raw/PREDICTIONS_R1.json'
    out = predict()
    if p.exists():
        old = json.load(open(p))
        assert out == old, 'Frozen prediction drift'
        print('R1 exact replay verified')
    else:
        write(p, out)
        now = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write(R / 'FROZEN_PREDICTIONS_R1.json', {'frozen_utc': now, 'prediction_file': 'raw/PREDICTIONS_R1.json', 'prediction_sha256': sha(p), 'prereg_sha256': sha(R / 'PREREG_R1.json'), 'existing_literature_numbers_seen': True, 'blinded_validation': False, 'future_local_measurement': False})
        print('R1 predictions frozen', len(out['rows']), now)
    write(R / 'raw/PREDICT_R1_COST.json', {'wall_seconds': time.perf_counter() - t, 'fit_seconds': 0})
