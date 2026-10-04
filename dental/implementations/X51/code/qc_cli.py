"""Freeze paired regional metrology prediction limits and score weekly observations."""
import thread_budget
import argparse, csv, datetime, json, math
from pathlib import Path
import numpy as np
from scipy.stats import t
from common import sha, clean
from drift import interval_hold

def dump(p, o):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(o), indent=2, allow_nan=False) + '\n')

def csvrows(path):
    return list(csv.DictReader(Path(path).open()))

def extract(rows, protocol):
    groups = {}
    for r in rows:
        if r.get('unit') != 'um' or not r.get('locator'):
            raise ValueError('Signed um values and locator required')
        if r.get('protocol_sha256') != protocol['protocol_sha256']:
            raise ValueError('Protocol changed')
        if r.get('frame_id') != protocol['frame_id']:
            raise ValueError('Reference frame mismatch')
        g = r['region']
        ident = r['measurement_id']
        ref = float(r['reference_um'])
        coupon = float(r['coupon_um'])
        if not np.isfinite(ref) or not np.isfinite(coupon):
            raise ValueError('Nonfinite metrology')
        if g not in protocol['regions']:
            raise ValueError('Unplanned region')
        dest = groups.setdefault(g, {})
        if ident in dest:
            raise ValueError('Duplicate measurement in region')
        dest[ident] = (ref, coupon - ref)
    if set(groups) != set(protocol['regions']):
        raise ValueError('Missing region')
    return {g: np.array(list(q.values())) for (g, q) in groups.items()}

def freeze(baseline, protocol, output):
    required = ['protocol_sha256', 'frame_id', 'regions', 'weeks', 'weekly_n', 'alpha', 'source_kind', 'slope_bound_um_per_day', 'hard_error_d_um', 'gap0_um', 'gap_limit_um', 'calibration_locator', 'model_scope']
    if any((k not in protocol for k in required)):
        raise ValueError('Incomplete QC protocol')
    if protocol['model_scope'] != 'COMMON_ADDITIVE_SIGNED_REGION':
        raise ValueError('Only declared additive signed regional model supported')
    if len(protocol['protocol_sha256']) != 64 or any((c not in '0123456789abcdef' for c in protocol['protocol_sha256'])):
        raise ValueError('Protocol SHA256 required')
    if len(set(protocol['regions'])) != len(protocol['regions']) or not protocol['regions']:
        raise ValueError('Unique region IDs required')
    if not isinstance(protocol['weeks'], int) or not 1 <= protocol['weeks'] <= 52 or (not isinstance(protocol['weekly_n'], int)) or (protocol['weekly_n'] < 2) or (not 0 < protocol['alpha'] < 1):
        raise ValueError('Invalid monitoring contract')
    if Path(output).exists():
        raise ValueError('Frozen prediction exists; start separately budgeted plan')
    groups = extract(baseline, protocol)
    R = len(groups)
    out = {}
    for (g, x) in groups.items():
        N = len(x)
        if N < 10:
            raise ValueError('At least10 independent baseline pairs per region required for this lab template')
        sd = x.std(axis=0, ddof=1)
        if np.any(sd <= 0):
            raise ValueError('Degenerate baseline cannot estimate variance')
        tail = protocol['alpha'] / (4 * R * protocol['weeks'])
        crit = float(np.nextafter(t.isf(tail, N - 1), np.inf))
        while t.sf(crit, N - 1) > tail:
            crit = float(np.nextafter(crit, np.inf))
        h = crit * sd * np.sqrt(1 / protocol['weekly_n'] + 1 / N)
        out[g] = {'N': N, 'mean_b_d_um': x.mean(axis=0).tolist(), 'sd_b_d_um': sd.tolist(), 'tcrit': crit, 'delta_limits_um': h.tolist()}
    o = {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'protocol': protocol, 'baseline': out, 'statistical_scope': 'Normal iid stable paired marginal observations; Bonferroni fixed regions x52 checks. No chart-driven reset/recalibration or extra looks under same budget.', 'physical_measurement_claim': protocol['source_kind'] == 'LAB_MEASUREMENT', 'standard_conformity': 'UNKNOWN'}
    dump(output, o)
    Path(str(output) + '.sha256').write_text(sha(output) + '\n')
    return o

def score(predfile, records, week, statefile):
    if sha(predfile) != Path(str(predfile) + '.sha256').read_text().strip():
        raise ValueError('Frozen metrology prediction changed')
    pred = json.loads(Path(predfile).read_text())
    proto = pred['protocol']
    groups = extract(records, proto)
    if not isinstance(week, int) or week < 1 or week > proto['weeks']:
        raise ValueError('Week outside fixed annual plan')
    digest = __import__('hashlib').sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()
    statepath = Path(statefile)
    history = json.loads(statepath.read_text()) if statepath.exists() else {'prediction_sha256': sha(predfile), 'weeks': {}}
    if history['prediction_sha256'] != sha(predfile):
        raise ValueError('QC state belongs to another plan')
    if str(week) in history['weeks']:
        if history['weeks'][str(week)]['input_sha256'] != digest:
            raise ValueError('No optional replacement/repeated look within a week')
        return history['weeks'][str(week)]['output']
    if history['weeks'] and week != max((int(w) for w in history['weeks'])) + 1:
        raise ValueError('Chronological complete weekly plan required')
    if not history['weeks'] and week != 1:
        raise ValueError('Start with week1')
    regions = {}
    any_alarm = False
    any_hold = False
    for (g, x) in groups.items():
        if len(x) != proto['weekly_n']:
            raise ValueError('Wrong weekly replicate count')
        mean = x.mean(axis=0)
        base = np.array(pred['baseline'][g]['mean_b_d_um'])
        delta = mean - base
        hits = np.abs(delta) > np.array(pred['baseline'][g]['delta_limits_um'])
        any_alarm = any_alarm or bool(hits.any())
        if proto['source_kind'] != 'LAB_MEASUREMENT' or not proto['calibration_locator'] or proto['slope_bound_um_per_day'] is None or (proto['hard_error_d_um'] is None):
            guard = {'verdict': 'HOLD_UNKNOWN', 'reason': 'independent physical bound/transfer calibration missing, or own fixture'}
        else:
            guard = interval_hold(proto['gap0_um'][g], float(delta[1]), proto['hard_error_d_um'], proto['slope_bound_um_per_day'], 7, proto['gap_limit_um'])
        any_hold = any_hold or guard['verdict'] != 'MODEL_WITHIN_LIMIT'
        regions[g] = {'delta_b_d_um': delta.tolist(), 'scanner_alarm': bool(hits[0]), 'process_alarm': bool(hits[1]), 'horizon_guard': guard, 'resolution': 'PER_SURFACE_REGION'}
    verdict = 'HOLD_CHART_ALARM' if any_alarm else 'HOLD_UNKNOWN_OR_HORIZON' if any_hold else 'CONDITIONAL_MODEL_WITHIN_LIMIT'
    o = {'week': week, 'verdict': verdict, 'regions': regions, 'prediction_sha256': sha(predfile), 'physical_release_authority': 'LAB_REQUIRED; verify evidence bounds and common additive model; no clinical or standards certification', 'next_operation': 'Check scanner against independent reference and region-preserving coupon/seated gap metrology; freeze any new model before new measurements'}
    history['weeks'][str(week)] = {'input_sha256': digest, 'output': o}
    dump(statepath, history)
    return o

def main():
    p = argparse.ArgumentParser(description=__doc__)
    s = p.add_subparsers(dest='cmd', required=True)
    q = s.add_parser('freeze')
    q.add_argument('--baseline', required=True)
    q.add_argument('--protocol', required=True)
    q.add_argument('--output', required=True)
    q = s.add_parser('score')
    q.add_argument('--prediction', required=True)
    q.add_argument('--measurements', required=True)
    q.add_argument('--week', type=int, required=True)
    q.add_argument('--state', required=True)
    a = p.parse_args()
    if a.cmd == 'freeze':
        o = freeze(csvrows(a.baseline), json.loads(Path(a.protocol).read_text()), a.output)
    else:
        o = score(a.prediction, csvrows(a.measurements), a.week, a.state)
    print(json.dumps(clean(o), indent=2, allow_nan=False))
if __name__ == '__main__':
    main()
