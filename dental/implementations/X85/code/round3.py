from common import *
from regional import predict
from geometry import roof_checks, ball_access, shell, export_stl
from round2 import load_item, simulated_calibration
from scipy.optimize import linprog
import time

def joint(model, B, weights, allow):
    J = np.array(model['J_N_per_mm'])
    f0 = np.array(model['f0'])
    beta = np.array(model['beta_N_per_mm'])
    eta = model['eta_N']
    G = []
    b = []
    for i in range(2):
        G.extend([np.r_[-J[i] + beta, 0.0], np.r_[J[i] + beta, 0.0]])
        b.extend([9 - eta - f0[i], f0[i] - eta - 7])
    G.extend(np.c_[J + beta[None, :], np.zeros(6)])
    b.extend(f0 - eta - 1e-06)
    G.extend(np.c_[B, -np.ones(len(B))])
    b.extend(allow)
    G = np.asarray(G)
    b = np.asarray(b)
    bounds = [(0, 0.05), (0, 0.05), (0, 2.0)]
    c = np.array([0.0, 0.0, 1.0])
    sol = linprog(c, A_ub=G, b_ub=b, bounds=bounds, method='highs', options={'threads': 1})
    if not sol.success:
        return dict(status='INFEASIBLE_OR_UNKNOWN', reason=sol.message)
    d = float(sol.x[2])
    dual = sol.ineqlin.marginals
    dl = sol.lower.marginals
    du = sol.upper.marginals
    objective_dual = float(b @ dual + np.array([0.0, 0.0, 0.0]) @ dl + np.array([0.05, 0.05, 2.0]) @ du)
    residual = float(np.max(abs(c - G.T @ dual - dl - du)))
    gap = float(abs(sol.fun - objective_dual))
    second = linprog(np.r_[np.asarray(weights) @ B, 0.0], A_ub=G, b_ub=b, bounds=[(0, 0.05), (0, 0.05), (d, d)], method='highs', options={'threads': 1})
    if not second.success:
        return dict(status='UNKNOWN_SECOND_STAGE', reason=second.message)
    below = None
    if d > 0.001:
        probe = linprog(np.zeros(3), A_ub=G, b_ub=b, bounds=[(0, 0.05), (0, 0.05), (0, d - 0.001)], method='highs', options={'threads': 1})
        below = dict(preparation_cap_mm=d - 0.001, solver_status=int(probe.status), infeasible=probe.status == 2)
    a = -second.x[:2]
    answer = predict(model, a)
    return dict(status='CONDITIONAL_MODEL_PREPARATION' if residual <= 1e-07 and gap <= 1e-07 else 'UNKNOWN_DUAL', minimum_extra_preparation_mm=d, baseline_wall_repair_mm=float(max(0, -np.min(allow))), increment_over_baseline_wall_repair_mm=d - max(0, -np.min(allow)), coefficients_mm=a, at_selected_height=answer, dual_stationarity_residual=residual, primal_dual_objective_gap_mm=gap, below_minimum_probe=below, geometry_changed_calibration='INVALIDATED_FOR_PHYSICAL_USE; fixed support retained only in simulation')

def run():
    start = time.perf_counter()
    rows = []
    q = v6_contact()
    for r0 in read(ROOT / 'raw/R2_ROWS.json'):
        if r0['status'] != 'EVALUATED':
            continue
        (t, z, inner) = load_item(r0)
        B = np.array(np.load(r0['artifact_path'])['basis'])
        (model, probes, A, C, P, f0) = simulated_calibration(t, z, inner, B)
        out = joint(model, B, t['weights'], z - inner - t['requirements']['wall_mm'])
        r = dict(uid=r0['uid'], joint=out)
        if out['status'] == 'CONDITIONAL_MODEL_PREPARATION':
            final = z + B @ out['coefficients_mm']
            inn = inner - out['minimum_extra_preparation_mm']
            r['geometry'] = roof_checks(t, final, inn)
            r['milling'] = ball_access(t['xy'], t['faces'], final)
            dest = ROOT / 'exports' / ('R3_' + r0['uid'])
            dest.mkdir(parents=True, exist_ok=True)
            (v, f) = shell(t['xy'], t['faces'], final, inn)
            export_stl(dest / 'crown_roof.stl', v, f)
            path = DATA / ('R3_' + r0['uid'] + '.npz')
            np.savez_compressed(path, xy=t['xy'], faces=t['faces'], outer_z=final, inner_z=inn, preparation_z=t['preparation_z'] - out['minimum_extra_preparation_mm'])
            r['artifact_path'] = str(path)
            r['artifact_sha256'] = sha(path)
            with np.load(V4 / 'payload/private/references' / f"{r0['case_key']}.npz", allow_pickle=False) as ref:
                rg = t['ceiling'] - ref[r0['family']]
            r['contact'] = q.compare(t['xy'], t['faces'], t['ceiling'] - final, rg)
        rows.append(r)
    summary = dict(requested=len(rows), joint_wall_force_candidates=sum((r['joint']['status'] == 'CONDITIONAL_MODEL_PREPARATION' and r['geometry']['wall_pass'] for r in rows)), all_geometry_passes=sum((r.get('geometry', {}).get('wall_pass', False) and r.get('geometry', {}).get('antagonist_pass', False) and (r.get('milling', {}).get('status') == 'CONDITIONAL_IDEAL_TOOL_PASS') for r in rows)), seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R3_ROWS.json', rows)
    dump(ROOT / 'raw/R3_SUMMARY.json', summary)
    state('R3_DECIDED', str(summary), 'Change surface representation to an ideal-tool accessible plane; measure source contact loss')
    (ROOT / 'history/HANDOFF_R3.md').write_text(json.dumps(summary) + '\nThe new die invalidates physical calibration; only a fixed-support scenario is solved. Minimum extra preparation is a model tradeoff, never clinical advice. Next: change actual upper-surface representation so the ideal milling proof can pass, retain external contact test.\n')
    print(json.dumps(summary))
if __name__ == '__main__':
    run()
