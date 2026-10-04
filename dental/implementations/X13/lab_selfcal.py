"""Freeze a prospective local three-region pilot from two calibration settings.

Inputs are LOCAL lab CSVs. This script never reads evaluation measurements.
The generated predictions require external held-dose validation before use.
"""
from pathlib import Path
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '4'
import csv, json, hashlib, datetime, argparse
import numpy as np

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def freeze(calibration, design, out):
    c = list(csv.DictReader(Path(calibration).open()))
    designrows = list(csv.DictReader(Path(design).open()))
    pred = []
    required = {'system', 'material', 'state', 'region', 'spacer_um', 'mean_um', 'n_specimens'}
    assert c and required <= set(c[0])
    assert designrows and (not {'mean_um', 'measured_mean_um'} & set(designrows[0])), 'Evaluation design must contain NO measured outcomes'
    for t in designrows:
        matched = [d for d in c if all((d[k] == t[k] for k in ['system', 'material', 'state', 'region']))]
        assert len(matched) == 2
        s = np.array([float(d['spacer_um']) for d in matched])
        y = np.array([float(d['mean_um']) for d in matched])
        assert len(set(s)) == 2 and np.all(s > 0)
        reciprocal = t['region'] == 'marginal'
        X = np.stack([np.ones(2), 1 / s if reciprocal else s], axis=1)
        b = np.linalg.solve(X, y)
        st = float(t['spacer_um'])
        mu = float(b[0] + b[1] * (1 / st if reciprocal else st))
        pred.append(dict(**t, predicted_group_mean_um=mu, basis='reciprocal' if reciprocal else 'affine', extrapolation=not min(s) <= st <= max(s), status='UNVALIDATED_LOCAL_PILOT', predictive_interval_um=None, physical_seating_prediction=None))
    f = Path(out)
    assert not f.exists(), 'Never overwrite a prospective prediction freeze'
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), calibration_path=str(Path(calibration).resolve()), calibration_sha256=sha(calibration), evaluation_design_path=str(Path(design).resolve()), evaluation_design_sha256=sha(design), predictions=pred, review_state='PENDING_INDEPENDENT_REVIEW', clinical_recommendation=False)
    f.write_text(json.dumps(record, indent=2) + '\n')
    f.with_suffix('.sha256').write_text(sha(f) + '\n')
if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('--calibration', required=True)
    a.add_argument('--test-design', required=True)
    a.add_argument('--output', default='FROZEN_LAB_PREDICTIONS.json')
    x = a.parse_args()
    freeze(x.calibration, x.test_design, x.output)
