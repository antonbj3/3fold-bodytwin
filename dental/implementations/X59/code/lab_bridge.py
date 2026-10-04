"""Actual-lab four-observation API. Does not manufacture calibration measurements.

freeze --calibration actual_pairs.csv --design future_replica_only.csv --output frozen.json
score --predictions frozen.json --validation later_actual_pairs.csv --output validation.json
"""
import argparse, csv, hashlib, json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
OBS = ['dry_ct_um', 'in_situ_replica_ct_um', 'sectioned_replica_ct_um', 'sectioned_replica_optical_um']
KEYS = ['specimen_id', 'protocol_id', 'region', 'point_id']

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load_csv(p):
    with Path(p).open() as f:
        return list(csv.DictReader(f))

def validate_rows(rows):
    if not rows:
        raise ValueError('No actual paired calibration measurements')
    seen = set()
    for row in rows:
        key = tuple((row.get(k, '') for k in KEYS))
        if not all(key) or key in seen:
            raise ValueError('Missing or duplicate specimen/protocol/region/point identity')
        seen.add(key)
        if not row.get('batch_id', ''):
            raise ValueError('Manufacturing batch identity is required separately from protocol')
        if row.get('point_pair_verified', '').lower() not in ['true', '1']:
            raise ValueError('Actual homologous point registration must be verified')
        if row.get('ct_state_bias_calibrated', '').lower() not in ['true', '1']:
            raise ValueError('Common CT state bias closure is unverified: require state-specific phantom calibration')
        for col in OBS:
            if col not in row:
                raise ValueError('Missing actual observation: ' + col)
            x = float(row[col])
            if not np.isfinite(x) or x < 0:
                raise ValueError('Invalid nonnegative gap measurement')
        bound = float(row.get('single_observation_bound_um', 'nan'))
        if not np.isfinite(bound) or bound < 0:
            raise ValueError('Missing measured observation/registration error bound')

def aggregate(rows):
    validate_rows(rows)
    groups = defaultdict(list)
    for r in rows:
        groups[tuple((r[k] for k in KEYS[:3]))].append(r)
    out = []
    for (key, rr) in groups.items():
        batches = {r['batch_id'] for r in rr}
        if len(batches) != 1:
            raise ValueError('One specimen/region cannot silently mix manufacturing batches')
        d = dict(zip(KEYS[:3], key))
        d.update({c: float(np.mean([float(r[c]) for r in rr])) for c in OBS})
        d['batch_id'] = next(iter(batches))
        d['region_difference_error_bound_um'] = 2 * float(np.mean([float(r['single_observation_bound_um']) for r in rr]))
        d['n_points'] = len(rr)
        d['resolution'] = 'PER_SURFACE_REGION'
        out.append(d)
    return out

def fit(rows):
    regional = aggregate(rows)
    groups = defaultdict(list)
    for r in regional:
        groups[r['protocol_id'], r['region']].append(r)
    out = []
    for ((protocol, region), rr) in groups.items():
        n = len(rr)
        if n < 10:
            raise ValueError('At least ten actual independent copings per protocol/region, not ten points')
        dif = np.array([r['sectioned_replica_optical_um'] - r['dry_ct_um'] for r in rr])
        bias = float(dif.mean())
        sd = float(dif.std(ddof=1))
        loo_resid = dif - (dif.sum() - dif) / (n - 1)
        loo_mae = float(np.mean(np.abs(loo_resid)))
        transitions = {name: float(np.mean([r[a] - r[b] for r in rr])) for (name, a, b) in [('seating', 'in_situ_replica_ct_um', 'dry_ct_um'), ('handling', 'sectioned_replica_ct_um', 'in_situ_replica_ct_um'), ('detector', 'sectioned_replica_optical_um', 'sectioned_replica_ct_um')]}
        out.append(dict(protocol_id=protocol, region=region, n_copings=n, calibration_batches=sorted({r['batch_id'] for r in rr}), bias_replica_minus_dryct_um=bias, sample_sd_paired_difference_um=sd, nominal_normal_loa_um=[bias - 1.96 * sd, bias + 1.96 * sd], normality_closure='Unverified; these are nominal sample LoA, not rigorous population coverage', leave_one_coping_out_mae_um=loo_mae, calibration_20um_gate='PASS' if loo_mae <= 20 else 'FAIL', source_mean_difference_bound_um=float(np.mean([r['region_difference_error_bound_um'] for r in rr])), isolated_contrasts_um=transitions, resolution='PER_SURFACE_REGION', empirical_new_batch_prediction_bound_um=None, absolute_true_gap_um=None, independent_specimen_assumption='User must ensure independently manufactured copings; IDs alone do not establish independence'))
    return out

def freeze(calibration, design, output):
    outpath = Path(output)
    if outpath.exists():
        raise ValueError('Prediction freeze already exists; preserve it')
    calibration_rows = load_csv(calibration)
    models = fit(calibration_rows)
    design_rows = load_csv(design)
    if not design_rows:
        raise ValueError('No future design/replica observations')
    train_ids = {r['specimen_id'] for r in calibration_rows}
    pred = []
    seen = set()
    for r in design_rows:
        if not r.get('batch_id', ''):
            raise ValueError('Future batch identity is required')
        if any(('ct' in c.lower() or 'validation' in c.lower() or 'observed_truth' in c.lower() for c in r)):
            raise ValueError('Future design CSV must not contain CT/validation outcomes')
        key = (r['specimen_id'], r['protocol_id'], r['region'])
        if key in seen or key[0] in train_ids:
            raise ValueError('Duplicate or training specimen in holdout')
        seen.add(key)
        model = next((m for m in models if m['protocol_id'] == key[1] and m['region'] == key[2]), None)
        if model is None or model['calibration_20um_gate'] != 'PASS':
            raise ValueError('No passing same-protocol regional calibration')
        v = float(r['replica_mean_um'])
        if not np.isfinite(v) or v < 0:
            raise ValueError('Invalid observed future replica mean')
        pred.append(dict(specimen_id=key[0], protocol_id=key[1], region=key[2], batch_id=r['batch_id'], replica_mean_um=v, predicted_dry_ct_mean_um=v - model['bias_replica_minus_dryct_um'], resolution='PER_SURFACE_REGION', prediction_bound_um=None, status='CONDITIONAL_METHOD_SCALE_PREDICTION_REQUIRES_HELDOUT_VALIDATION'))
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), claim_type='information_link', calibration_path=str(Path(calibration).resolve()), calibration_sha256=sha(calibration), design_path=str(Path(design).resolve()), design_sha256=sha(design), models=models, predictions=pred, disclosure='Actual input CSVs, pending future independent validation. CT is a method comparator; absolute gap remains unknown.')
    outpath.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    outpath.with_suffix(outpath.suffix + '.sha256').write_text(sha(outpath) + '\n')
    return record

def score(predictions, validation, output):
    path = Path(predictions)
    if sha(path) != path.with_suffix(path.suffix + '.sha256').read_text().strip():
        raise ValueError('Frozen predictions changed')
    f = json.loads(path.read_text())
    actual = aggregate(load_csv(validation))
    scores = []
    for p in f['predictions']:
        a = next((r for r in actual if all((r[k] == p[k] for k in KEYS[:3]))), None)
        if a is None:
            raise ValueError('Missing held-out validation specimen/region')
        if a['batch_id'] != p['batch_id']:
            raise ValueError('Validation batch differs from frozen design identity')
        if abs(a['sectioned_replica_optical_um'] - p['replica_mean_um']) > a['region_difference_error_bound_um']:
            raise ValueError('Future replica mean and validation optical observation state disagree beyond supplied bound')
        error = abs(a['dry_ct_um'] - p['predicted_dry_ct_mean_um'])
        scores.append(dict(**p, measured_dry_ct_mean_um=a['dry_ct_um'], absolute_error_um=error, observed_error_upper_bound_um=error + a['region_difference_error_bound_um'], gate_20um='PASS' if error + a['region_difference_error_bound_um'] <= 20 else 'FAIL', comparator_raw_replica_error_um=abs(a['dry_ct_um'] - p['replica_mean_um'])))
    result = dict(predictions_sha256=sha(path), validation_path=str(Path(validation).resolve()), validation_sha256=sha(validation), scores=scores, overall_gate='PASS' if all((s['gate_20um'] == 'PASS' for s in scores)) else 'FAIL', status='Target-protocol held-out observation comparison; future batches and absolute true gap remain unvalidated')
    if Path(output).exists():
        raise ValueError('Preserve previous validation version')
    Path(output).write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('freeze')
    a.add_argument('--calibration', required=True)
    a.add_argument('--design', required=True)
    a.add_argument('--output', required=True)
    b = sub.add_parser('score')
    b.add_argument('--predictions', required=True)
    b.add_argument('--validation', required=True)
    b.add_argument('--output', required=True)
    args = vars(parser.parse_args())
    cmd = args.pop('cmd')
    result = freeze(**args) if cmd == 'freeze' else score(**args)
    print(json.dumps(dict(status=result.get('overall_gate', 'FROZEN'), output=args['output'])))
if __name__ == '__main__':
    main()
