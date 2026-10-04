"""An explicit abstaining inverse API; no clinical recommendations."""
from pathlib import Path
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '4'
import json, time, argparse
import numpy as np
P = Path(__file__).resolve().parent

def regional_query(region, low, high):
    r = json.loads((P / 'RESULTS_R1.json').read_text())['metrics'][region]
    return {'region': region, 'target_um': [low, high], 'CAD_setting_um': None, 'status': 'UNKNOWN', 'reason': 'Insufficient independently measured within-system spacer perturbations; predictive error/coverage gate failed', 'support': r['dose_support'], 'model_error_um': r['candidate_rmse_um'], 'clinical_recommendation': False}

def local_marginal_query(study, low, high):
    m = next((d for d in json.loads((P / 'FROZEN_PREDICTIONS_R3.json').read_text())['models'] if d['study'] == study))
    s = np.linspace(*m['interpolation_domain_um'], 1001)
    mu = m['b_um'] + m['a_um_squared'] / s
    ok = (mu >= low) & (mu <= high)
    return dict(study=study, region='marginal', target_um=[low, high], status='CONDITIONAL_LOCAL_PILOT', central_curve_in_target_CAD_range_um=[float(s[ok].min()), float(s[ok].max())] if ok.any() else None, calibrated_domain_CAD_um=m['interpolation_domain_um'], certified=False, uncertainty='Source anchor SD and covariance do not establish a measured predictive interval', extrapolation_allowed=False, clinical_recommendation=False, caution='Applies to the same study system/material and dry vertical marginal group mean; does not predict axial/occlusal film or a new batch')

def observed_decisions(rows, targets):
    arms = {}
    for d in rows:
        if d['region'] not in targets:
            continue
        key = (d['study'], d['arm'])
        arms.setdefault(key, {})[d['region']] = d
    out = []
    for ((study, arm), rs) in sorted(arms.items()):
        if set(rs) != set(targets):
            continue
        vals = {k: rs[k]['measured_mean_um'] for k in targets}
        passr = {k: targets[k][0] <= v <= targets[k][1] for (k, v) in vals.items()}
        out.append(dict(study=study, arm=arm, internal_spacer_um=rs['axial']['internal_spacer_um'], observed_mean_um=vals, within_engineering_target=passr, all_group_means_in_target=all(passr.values()), status='DESCRIPTIVE_OBSERVED_GROUP_MEANS_ONLY', n_specimens=rs['axial']['n_specimens']))
    return out

def demo():
    from predict_r1 import dump
    t = time.perf_counter()
    targets = json.loads((P / 'PREREG_R1.json').read_text())['illustrative_targets_um']
    rows = json.loads((P / 'measurements.json').read_text())
    out = dict(regional=[regional_query(k, *v) for (k, v) in targets.items()], local_marginal=[local_marginal_query(s, *targets['marginal']) for s in ['PMC10246932', 'PMC10721348']], observed_joint=observed_decisions(rows, targets), seconds=time.perf_counter() - t, scope='Engineering demo targets; no clinical advice')
    dump(P / 'INVERSE_DECISIONS.json', out)
    return out
if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('--region', choices=['marginal', 'axial', 'occlusal'], default='axial')
    a.add_argument('--target', nargs=2, type=float, default=[50, 100])
    a.add_argument('--local-study', choices=['PMC10246932', 'PMC10721348'])
    args = a.parse_args()
    print(json.dumps(local_marginal_query(args.local_study, *args.target) if args.local_study else regional_query(args.region, *args.target), indent=2))
