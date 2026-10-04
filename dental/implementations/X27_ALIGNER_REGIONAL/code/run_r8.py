import json, hashlib, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from whole_arch_refined import build, solve_floating
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R8.json').read_text())
    ge = pr['gates']
    assert hashlib.sha256((R / 'PREREG_R8.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R8.sha256.json').read_text())['sha256']
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    rows = []
    failures = []
    for arch in arches:
        s = build(arch, regions)
        for active in [11, 13, 15]:
            if str(active) not in arch['teeth']:
                continue
            try:
                cand = solve_floating(s, active)
                ctrl = solve_floating(s, active, 'primal')
                parity = max((float(np.max(np.abs(np.array(cand['wrenches'][k]) - ctrl['wrenches'][k]))) for k in cand['wrenches']))
                rows.append(dict(case=arch['case'], active_fdi=active, source_teeth=len(s['ids']), candidate=cand, control=ctrl, max_wrench_difference=parity))
            except RuntimeError as ex:
                failures.append(dict(case=arch['case'], active_fdi=active, error=str(ex)))
    central = []
    for r in rows:
        if r['active_fdi'] != 11:
            continue
        a = next((a for a in arches if a['case'] == r['case']))
        e = np.array(a['teeth']['11']['buccal_unit'])
        central.append(float(np.array(r['candidate']['wrenches']['11'][:3]) @ e))
    gamma = 1.49 / float(np.median(central)) if central and np.median(central) > 1e-12 else None
    payload = dict(rows=rows, failures=failures, central_gain=gamma)
    h = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R8.json').read_bytes()).hexdigest(), predictions_sha256=h, predictions=payload)
    if (R / 'FROZEN_PREDICTIONS_R8.json').exists():
        if json.loads((R / 'FROZEN_PREDICTIONS_R8.json').read_text())['predictions_sha256'] != h:
            raise RuntimeError('R8 drift')
    else:
        put('FROZEN_PREDICTIONS_R8.json', frozen)
    source = json.loads((R / 'inputs/THESIS_EXTERNAL_ADDENDUM.json').read_text())
    comparisons = []
    for ref in source['neighbours']:
        values = []
        for row in rows:
            if row['active_fdi'] != ref['active_fdi'] or gamma is None or str(ref['response_fdi']) not in row['candidate']['wrenches']:
                continue
            arch = next((a for a in arches if a['case'] == row['case']))
            e = np.array(arch['teeth'][str(ref['response_fdi'])]['buccal_unit'])
            values.append(float(gamma * np.array(row['candidate']['wrenches'][str(ref['response_fdi'])][:3]) @ e))
        prediction = float(np.median(values)) if values else None
        error = abs(prediction - ref['force_N']) / abs(ref['force_N']) if prediction is not None else None
        comparisons.append(dict(**ref, predicted_N=prediction, relative_error=error, gate=error is not None and error <= ge['neighbour_relative_error_max']))
    active_comp = []
    for (fdi, truth) in [(13, 2.25), (15, 1.5)]:
        vals = []
        for row in rows:
            if row['active_fdi'] != fdi or gamma is None:
                continue
            arch = next((a for a in arches if a['case'] == row['case']))
            e = np.array(arch['teeth'][str(fdi)]['buccal_unit'])
            vals.append(float(gamma * np.array(row['candidate']['wrenches'][str(fdi)][:3]) @ e))
        pred = float(np.median(vals)) if vals else None
        err = abs(pred - truth) / truth if pred is not None else None
        active_comp.append(dict(fdi=fdi, predicted_N=pred, external_N=truth, relative_error=err, gate=err is not None and err <= ge['unfit_active_relative_error_max']))
    gates = dict(neighbours=all((r['gate'] for r in comparisons)), unfit_active=all((r['gate'] for r in active_comp)), numerical_parity=bool(rows) and all((r['max_wrench_difference'] <= ge['primal_dual_wrench_abs'] for r in rows)), force_balance=bool(rows) and all((r['candidate']['force_balance_N'] <= ge['force_balance_N'] for r in rows)), torque_balance=bool(rows) and all((r['candidate']['torque_balance_Nmm'] <= ge['torque_balance_Nmm'] for r in rows)), nonpenetration=bool(rows) and all((r['candidate']['minimum_gap_mm'] >= -ge['nonpenetration_mm'] for r in rows)), null_modes=bool(rows) and all((r['candidate']['null_modes'] == 6 for r in rows)), injected_support_force_rejected=False)
    corrupted = []
    for row in rows:
        wrench = {k: np.array(v).copy() for (k, v) in row['candidate']['wrenches'].items()}
        wrench[str(row['active_fdi'])][0] += 0.01
        corrupted.append(float(np.linalg.norm(sum((v[:3] for v in wrench.values()), np.zeros(3)))) > ge['force_balance_N'])
    gates['injected_support_force_rejected'] = bool(corrupted) and all(corrupted)
    result = dict(round='R8', claim_type='algorithm', outcome='WHOLE_ARCH_CONTACT_FORCE_TRANSFER_PASS' if all(gates.values()) else 'WHOLE_ARCH_CONTACT_FORCE_TRANSFER_FAIL', external_referent=pr['external_referent'], gates=gates, central_gain=gamma, neighbour_rows=comparisons, active_rows=active_comp, failures=failures, prediction_file='FROZEN_PREDICTIONS_R8.json', wrench_scenarios=len(rows), wall_s=time.perf_counter() - start, maximum_numerical_wrench_difference=max((r['max_wrench_difference'] for r in rows)) if rows else None, next_operation='Same-shell measured geometry/contact/retained shape and nonlinear shell membrane law; whole-arch tooth-only support topology now replaces fictitious housing grounds')
    put('round8/results.json', result)
    (R / 'round8/HANDOFF.md').write_text('# R8 whole-arch tooth-supported shell\n\n' + result['outcome'] + '\n\n' + json.dumps(gates) + '\n\n' + result['next_operation'] + '\n')
    print(json.dumps(dict(outcome=result['outcome'], gates=gates, central_gain=gamma, neighbors=comparisons, active=active_comp, failures=failures), indent=2))
if __name__ == '__main__':
    main()
