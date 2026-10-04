import json, time, hashlib
import numpy as np
from scipy.optimize import linprog
from run_r1 import ROOT, write, freeze, state, gate

def distribute(x, upper, budget, order):
    x = x.copy()
    remaining = budget - x.sum()
    for j in order:
        add = min(max(0, upper[j] - x[j]), max(0, remaining))
        x[j] += add
        remaining -= add
    if abs(remaining) > 1e-07:
        raise RuntimeError('Cannot allocate budget')
    return x

def bounds(a, u, b):
    if b < 0 or b > u.sum() + 1e-08:
        raise ValueError('Published volume scenario not feasible in anterior-only box')
    lower_i = a + np.maximum(0, b - (u.sum() - u))
    idx = int(lower_i.argmin())
    xl = np.zeros(len(a))
    xl[idx] = lower_i[idx] - a[idx]
    xl = distribute(xl, u, b, [j for j in np.argsort(a)[::-1] if j != idx])
    lo = float(a.min())
    hi = float((a + u).min())
    for _ in range(70):
        t = (lo + hi) / 2
        if np.maximum(0, t - a).sum() <= b:
            lo = t
        else:
            hi = t
    xh = distribute(np.maximum(0, lo - a), u, b, np.argsort(a)[::-1])
    return (float((a + xl).min()), float((a + xh).min()), xl, xh)

def lp_control(a, u, b):
    n = len(a)
    q = np.zeros(n + 1)
    q[-1] = -1
    ub = np.c_[-np.eye(n), np.ones(n)]
    eq = np.r_[np.ones(n), 0.0][None, :]
    res = linprog(q, A_ub=ub, b_ub=a, A_eq=eq, b_eq=[b], bounds=[(0, float(x)) for x in u] + [(None, None)], method='highs')
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.x[-1])

def acquisition(a, u, b, upper, measurements, epsilon=10.0):
    ix = np.where(a < upper - epsilon)[0]
    observed = measurements[ix]
    other = np.ones(len(a), bool)
    other[ix] = False
    lower = min(float(observed.min()) if len(ix) else np.inf, float(a[other].min()) if other.any() else np.inf)
    high = min(upper, float(observed.min()) if len(ix) else np.inf)
    return (ix, float(lower), float(high))

def main():
    start = time.monotonic()
    p = ROOT / 'PREREG_R3.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / (p.name + '.sha256')).read_text().split()[0]
    ps = json.loads((ROOT / 'raw/profiles.json').read_text())
    cs = json.loads((ROOT / 'raw/cases_R1.json').read_text())
    refs = json.loads((ROOT / 'raw/round_R1.json').read_text())['external_comparisons']
    rows = []
    witnesses = {}
    max_res = max_err = max_acq = 0.0
    fast_s = lp_s = 0.0
    for c in cs:
        case = c['case']
        if case not in ps:
            continue
        pr = ps[case]
        elig = np.array(pr['eligible'], bool)
        a = np.array(pr['area_mm2'])[elig]
        w = np.array(pr['lateral_support_mm'])[elig]
        z = np.array(pr['z_mm'])[elig]
        dz = c['spacing_zyx_mm'][0]
        for ref in refs:
            t = time.monotonic()
            u = ref['dose_mm'] * w
            b = (ref['v1_cm3'] / ref['v0_cm3'] - 1) * a.sum()
            row = {'case': case, 'external_scenario_group': ref['group'], 'dose_mm': ref['dose_mm'], 'volume_ratio': ref['v1_cm3'] / ref['v0_cm3'], 'n_slices': len(a)}
            try:
                (lo, hi, xlo, xhi) = bounds(a, u, b)
            except ValueError as e:
                row.update(status='INFEASIBLE_CONDITIONAL_CLOSURE', reason=str(e))
                rows.append(row)
                continue
            fast_s += time.monotonic() - t
            t = time.monotonic()
            ctrl = lp_control(a, u, b)
            lp_s += time.monotonic() - t
            lower_i = a + np.maximum(0, b - (u.sum() - u))
            idx = int(lower_i.argmin())
            obj = np.zeros(len(a))
            obj[idx] = 1
            res = linprog(obj, A_eq=np.ones((1, len(a))), b_eq=[b], bounds=[(0, float(v)) for v in u], method='highs')
            if not res.success:
                raise RuntimeError(res.message)
            lowctrl = float(a[idx] + res.x[idx])
            residual = max(abs(xlo.sum() - b), abs(xhi.sum() - b)) * dz
            error = max(abs(ctrl - hi), abs(lo - lowctrl))
            (ix, l, h) = acquisition(a, u, b, hi, a + xlo)
            (_, l2, h2) = acquisition(a, u, b, hi, a + xhi)
            xu = b * u / u.sum()
            (_, l3, h3) = acquisition(a, u, b, hi, a + xu)
            aw = max(h - l, h2 - l2, h3 - l3)
            max_res = max(max_res, residual)
            max_err = max(max_err, error)
            max_acq = max(max_acq, aw)
            key = case + '_' + ref['group']
            witnesses[key] = {'area0_mm2': a.tolist(), 'area_low_mm2': (a + xlo).tolist(), 'area_high_mm2': (a + xhi).tolist(), 'z_mm': z.tolist(), 'upper_increment_mm2': u.tolist(), 'budget_area_sum_mm2': b, 'spacing_z_mm': dz, 'required_measurement_slice_indices': ix.tolist(), 'required_measurement_array_z_mm': z[ix].tolist()}
            row.update(status='FEASIBLE_CONDITIONAL', minimum_lower_mm2=lo, minimum_upper_mm2=hi, width_mm2=hi - lo, volume_witness_residual_mm3=residual, LP_difference_mm2=error, proportional_volume_point_mm2=float(a.min() * ref['v1_cm3'] / ref['v0_cm3']), required_slices=len(ix), required_slice_fraction=len(ix) / len(a), required_z_mm=z[ix].tolist(), acquisition_max_width_on_own_witnesses_mm2=aw, low_bottleneck_z_mm=float(z[np.argmin(a + xlo)]), high_bottleneck_z_mm=float(z[np.argmin(a + xhi)]))
            rows.append(row)
    frozen = freeze('FROZEN_PREDICTIONS_R3.json', {'intervals': rows, 'interpretation': 'Conditional feasible geometry extremes before any new measurement; no real post-TF2 measurement exists.'})
    good = [r for r in rows if r['status'] == 'FEASIBLE_CONDITIONAL']
    fraction = sum((gate(r['width_mm2'], 10.0) for r in good)) / len(good)
    med = float(np.median([r['required_slice_fraction'] for r in good]))
    controls = {'volume_witness': gate(max_res, 1e-06), 'LP_agreement': gate(max_err, 1e-06), 'scalar_volume_identifies': fraction >= 0.8, 'acquisition_guarantee': gate(max_acq, 10.0), 'compact_acquisition': gate(med, 0.25)}
    first = next(iter(witnesses.values()))
    mut = np.array(first['area_high_mm2'])
    mut[0] += 1.0
    wrong_volume = abs((mut - np.array(first['area0_mm2'])).sum() - first['budget_area_sum_mm2']) * first['spacing_z_mm']
    out = {'round': 'R3', 'outcome': 'SCALAR_VOLUME_REFUTED_WITH_CONDITIONAL_ACQUISITION' if not controls['scalar_volume_identifies'] else 'CONDITIONAL_VOLUME_IDENTIFIABLE', 'scenarios': rows, 'controls': controls, 'injected_error_rejected': {'volume_witness': not gate(wrong_volume, 1e-06), 'LP_agreement': not gate(max_err + 1.0, 1e-06), 'scalar_volume_identifies': 0.0 >= 0.8, 'acquisition_guarantee': not gate(20.0, 10.0), 'compact_acquisition': not gate(1.0, 0.25)}, 'feasible_scenarios': len(good), 'within_10mm2_fraction': fraction, 'median_required_slice_fraction': med, 'max_numeric_error_mm2': max_err, 'max_volume_residual_mm3': max_res, 'cost': {'execution_seconds': time.monotonic() - start, 'closed_form_bound_seconds': fast_s, 'LP_upper_seconds': lp_s, 'fit_parameters': 0, 'physical_measurements_performed': 0}, 'prediction_sha256': frozen, 'external_referent': json.loads(p.read_text())['external_referent'], 'acquisition_status': 'Proved sufficient set in the declared nonnegative expansion model; minimal-cardinality set not proved. Physical validity UNKNOWN.', 'facit_scope': 'TF2 expert labels are baseline measurement. Shi external dose/volume ratios define scenarios, not TF2 post-response targets. LP checks only the mathematics, not tissue mechanics.'}
    out['injected_error_rejected']['scalar_volume_identifies'] = not 0.0 >= 0.8
    write('raw/round_R3.json', out)
    write('raw/witnesses_R3.json', witnesses)
    state('R3 adjudicated', out['outcome'], 'Export runnable demo, figure, raw witnesses and scoped graph coverage proposal')
    print(json.dumps({k: v for (k, v) in out.items() if k != 'scenarios'}, indent=2))
if __name__ == '__main__':
    main()
