from pathlib import Path
import json, hashlib, datetime, time
from series import fit, stiffness
R = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(n, x):
    (R / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
x = json.load(open(R / 'inputs/R2_CALIBRATION_ENDPOINTS.json'))
rows = []
fits = {}
maxdiff = 0
t0 = time.perf_counter()
for (direction, ep) in x['calibration'].items():
    f = fit(ep)
    c = fit(ep, True)
    fits[direction] = f
    for s in [0.4, 0.5, 0.625]:
        k = stiffness(f, s)
        kk = stiffness(c, s)
        maxdiff = max(maxdiff, abs(k - kk) / k)
        rows.append({'direction': direction, 'sheet_mm': s, 'stiffness_N_per_mm': k, 'force_0.2mm_N': k * 0.2, 'force_0.25mm_N': k * 0.25, 'point_calibrated_cubic_K': ep[1]['K_N_per_mm'] * (s / 0.75) ** 3, 'constant_stiffness_K': ep[1]['K_N_per_mm'], 'primary': s == 0.625})
    fit1 = ep[1]['K_N_per_mm']
    fits[direction]['one_measurement_counterexample'] = {'same_measured_K_at_0.75': fit1, 'K_at_0.3_pure_bending': fit1 * (0.3 / 0.75) ** 3, 'K_at_0.3_saturated': fit1, 'identified_with_one_measurement': False}
out = {'schema': 'x19-predictions-r2-v1', 'claim_type': 'capability', 'fits': fits, 'rows': rows, 'control_max_relative_difference': maxdiff, 'prereg_sha256': sha(R / 'PREREG_R2.json'), 'calibration_sha256': sha(R / 'inputs/R2_CALIBRATION_ENDPOINTS.json'), 'code_sha256': {p.name: sha(p) for p in [R / 'code/series.py', R / 'code/predict_r2.py']}, 'scope': 'Within published sensorseries; cannot transport originalsheet parametersto measuredformedthickness or realarches without pairedcalibration'}
p = R / 'raw/PREDICTIONS_R2.json'
if p.exists():
    assert json.load(open(p)) == out, 'Frozen R2 prediction drift'
    print('R2 exact replay verified')
else:
    write('raw/PREDICTIONS_R2.json', out)
    write('FROZEN_PREDICTIONS_R2.json', {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prediction_file': 'raw/PREDICTIONS_R2.json', 'prediction_sha256': sha(p), 'primary_intermediate_validation_cell_not_read_yet': True, 'independent_lab_measurement': False, 'known_published_study': 'Yes; qualitative plateau used as prior', 'prereg_sha256': sha(R / 'PREREG_R2.json')})
    print(json.dumps(rows, indent=2))
write('raw/PREDICT_R2_COST.json', {'wall_seconds': time.perf_counter() - t0, 'endpoint_cells': 4, 'fits': 2})
