"""Prospective comparison against an unaltered same-pair measurement CSV.
This command does not create force observations or calibrate against the test.
"""
import argparse, csv, json, hashlib, datetime
from pathlib import Path
import numpy as np
from occlusion_operator import H, REGIONS, dump, sha

def main():
    a = argparse.ArgumentParser()
    a.add_argument('measurement_csv')
    a.add_argument('metadata_json')
    a.add_argument('--out', default='MEASURED_COMPARISON.json')
    args = a.parse_args()
    freeze = H / 'FROZEN_PREDICTIONS_R4.json'
    if not freeze.exists():
        raise RuntimeError('Predictions not frozen before acquisition')
    frozen = json.loads(freeze.read_text())
    meta = json.loads(Path(args.metadata_json).read_text())
    case = int(meta['case'])
    pred = H / 'raw/cases_3d' / f'{case:03d}.json'
    if sha(pred) != frozen['cases_sha256'][str(case)]:
        raise RuntimeError('Prediction hash drift')
    r = json.loads(pred.read_text())
    for field in ['upper_sha256', 'lower_sha256']:
        if meta[field] != r[field]:
            raise RuntimeError('Not the frozen registered geometry')
    if meta['quantity_unit'] not in ['N', 'pp']:
        raise ValueError('Unit must be N or pp')

    def utc_date(s):
        d = datetime.datetime.fromisoformat(s.replace('Z', '+00:00'))
        if d.tzinfo is None or d.utcoffset() != datetime.timedelta(0):
            raise ValueError('Timestamp must be explicit UTC')
        return d
    if utc_date(meta['acquired_utc']) <= utc_date(frozen['timestamp_utc']):
        raise ValueError('Measurement must postdate prediction freeze')
    for field in ['instrument', 'calibration_locator', 'observation_map']:
        if not meta.get(field):
            raise ValueError('Missing observation metadata: ' + field)
    observations = list(csv.DictReader(open(args.measurement_csv)))
    values = {}
    for row in observations:
        key = (row['jaw'], row['region'])
        v = float(row['value'])
        if row['jaw'] not in ['upper', 'lower'] or row['region'] not in REGIONS or key in values or (not np.isfinite(v)) or (v < 0):
            raise ValueError('Invalid/duplicate row')
        values[key] = v
    expected = {(j, r) for j in ['upper', 'lower'] for r in REGIONS}
    if set(values) != expected:
        raise ValueError('Need exact16paired regional channels; use documented observation map for sensor coordinates')
    out = []
    for (jaw, offset) in [('upper', 0), ('lower', 8)]:
        observed = np.array([values[jaw, region] for region in REGIONS])
        if observed.sum() <= 0:
            raise ValueError('Zero force/signal arch')
        if meta['quantity_unit'] == 'N':
            observed = observed / observed.sum() * 100
        elif abs(observed.sum() - 100) > 5:
            raise ValueError('Relative signal normalization inconsistent')
        metric = r['metrics'][1]
        w = np.asarray(metric['mechanics']['shares_pp'])[offset:offset + 8]
        area = np.asarray(metric['area_shares_pp'])[offset:offset + 8]
        out.append(dict(jaw=jaw, observed_pp=observed.tolist(), mechanics_MAE_pp=float(np.abs(observed - w).mean()), area_MAE_pp=float(np.abs(observed - area).mean()), gate='SAME_PATIENT_TEST_ONLY; anatomical mapping/sensor-transfer metadata require independent review'))
    dump(H / args.out, dict(frozen_sha256=sha(freeze), measurement_sha256=sha(args.measurement_csv), metadata_sha256=sha(args.metadata_json), rows=out, calibration_performed=False, review_state='PENDING_INDEPENDENT_REVIEW'))
if __name__ == '__main__':
    main()
