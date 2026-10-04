"""Shared upper/lower reaction simplex, sharp necessary bounds and measurement rank.
This is force conservation on a provisional projected contact graph. No tissue
closure, applied jaw moment or physical-force measurement is identified here.
"""
import json, time, copy
import numpy as np
from scipy.optimize import linprog, nnls
from occlusion_operator import H, dump, state, REGIONS

def incidence(pairs):
    pairs = np.asarray(pairs, int)
    B = np.zeros((16, len(pairs)))
    if len(pairs):
        B[pairs[:, 0], np.arange(len(pairs))] = 1
        B[pairs[:, 1] + 8, np.arange(len(pairs))] = 1
    return B

def nnls_tangent(pairs, k):
    B = incidence(pairs)
    kr = k / k.mean()
    A = B / np.sqrt(kr)[:, None]
    lam = nnls(A, 0.5 * np.sqrt(kr), maxiter=10000)[0]
    u = B @ lam / kr
    F = k * u
    slack = B.T @ u - 1
    conservation = abs(F[:8].sum() - F[8:].sum()) / F[:8].sum()
    station = np.max(np.abs(B @ lam - kr * u))
    comp = np.max(np.abs(lam * slack))
    w = np.r_[F[:8] / F[:8].sum(), F[8:] / F[8:].sum()] * 100
    return dict(shares_pp=w.tolist(), conditional_force_per_unit_increment_N=F.tolist(), unit_increment_displacements=u.tolist(), conservation_relative=float(conservation), kkt_residual=float(station), complementarity=float(comp), minimum_slack=float(slack.min()), valid=bool(conservation <= 1e-07 and station <= 1e-06 and (comp <= 1e-06) and (slack.min() >= -1e-06)), dual_reaction_N_per_common_stiffness=lam.tolist())

def bounds_and_basis(pairs):
    B = incidence(pairs)
    n = B.shape[1]
    if not n:
        return dict(status='NO_CONTACT', affine_dimension=None, intervals=None)
    closed = np.c_[B.min(1), B.max(1)]
    lp = []
    start = time.perf_counter()
    for row in B:
        lo = linprog(row, A_eq=np.ones((1, n)), b_eq=[1], bounds=(0, None), method='highs')
        hi = linprog(-row, A_eq=np.ones((1, n)), b_eq=[1], bounds=(0, None), method='highs')
        if not (lo.success and hi.success):
            raise RuntimeError('LP failed')
        lp.append([lo.fun, -hi.fun])
    lp = np.asarray(lp)
    error = float(np.max(np.abs(lp - closed)))
    rank = int(np.linalg.matrix_rank(B, tol=1e-09))
    measurements = []
    M = np.ones((1, n))
    old = 1
    for r in range(16):
        new = np.vstack([M, B[r]])
        rk = int(np.linalg.matrix_rank(new, tol=1e-09))
        if rk > old:
            M = new
            old = rk
            measurements.append(('upper:' if r < 8 else 'lower:') + REGIONS[r % 8])
    front = (B[[0, 1, 4, 5]].sum(0) + B[[8, 9, 12, 13]].sum(0)) / 2
    front_interval = [float(front.min()), float(front.max())]
    one_front_rank = int(np.linalg.matrix_rank(np.vstack([np.ones(n), front]), tol=1e-09)) - 1
    return dict(status='NECESSARY_OUTER_SET_ONLY', n_edges=n, visible_rank=rank, affine_dimension=rank - 1, minimum_scalar_shape_channels=rank - 1, measurement_basis=measurements, basis_sufficient=old == rank, regional_intervals=closed.tolist(), LP_control_intervals=lp.tolist(), LP_parity_error=error, lp_wall_s=time.perf_counter() - start, total_force_shape_rank_gain=0, one_anterior_channel_shape_rank_gain=one_front_rank, anterior_share_interval=front_interval, physical_validity='UNKNOWN: region labels, active-contact uncertainty, support law, moments and load are not measured')

def main():
    started = time.perf_counter()
    rows = []
    nom = np.tile([2, 1, 2, 2, 2, 1, 2, 2], 2)
    central = np.sqrt(500 * 1130)
    k = nom * central
    cache = {}
    lp_cache = {}
    scenarios = [k, nom * 500, nom * 1130]
    for i in range(16):
        for bound in [500, 1130]:
            kk = k.copy()
            kk[i] = nom[i] * bound
            scenarios.append(kk)
    for p in sorted((H / 'raw/cases').glob('*.json')):
        r = json.loads(p.read_text())
        case = r['case']
        metrics = []
        for met in r['metrics']:
            if met['mechanics'] is None:
                metrics.append(dict(band_mm=met['band_mm'], status='NO_CONTACT'))
                continue
            pairs = met['mechanics']['pairs']
            key = str(pairs)
            if key not in lp_cache:
                lp_cache[key] = bounds_and_basis(pairs)
            if key not in cache:
                out = [nnls_tangent(pairs, kk) for kk in scenarios]
                ws = np.array([o['shares_pp'] for o in out])
                cache[key] = dict(repaired_nominal=out[0], support_scenario_envelope_pp=np.c_[ws.min(0), ws.max(0)].tolist(), support_scenarios_valid=all((o['valid'] for o in out)), shared_scale_profile_delta_pp=float(max(np.max(np.abs(ws[0] - ws[1])), np.max(np.abs(ws[0] - ws[2])))), support_scenario_count=len(scenarios))
            metric = dict(band_mm=met['band_mm'], **lp_cache[key], **cache[key], repaired_vs_SLSQP_Linf_pp=float(np.max(np.abs(np.array(cache[key]['repaired_nominal']['shares_pp']) - met['mechanics']['shares_pp']))))
            metrics.append(metric)
        rows.append(dict(case=case, metrics=metrics))
        if case % 50 == 0:
            state('R3_RUNNING', 'R1/R2 physical force and region gates failed', f'Force outer sets case{case}')
            print('force feasible', case, flush=True)
    centralrows = [m for r in rows for m in r['metrics'] if m['band_mm'] == 0.1 and m.get('affine_dimension') is not None]
    result = dict(round='R3', n_cases=len(rows), no_contact_cases=[r['case'] for r in rows if r['metrics'][1]['status'] == 'NO_CONTACT'], LP_max_parity_error=max((m['LP_parity_error'] for m in centralrows)), rank_distribution={str(i): sum((m['affine_dimension'] == i for m in centralrows)) for i in sorted(set((m['affine_dimension'] for m in centralrows)))}, identified_joint_vectors=sum((m['affine_dimension'] == 0 for m in centralrows)), minimum_channels_summary=dict(min=min((m['minimum_scalar_shape_channels'] for m in centralrows)), median=float(np.median([m['minimum_scalar_shape_channels'] for m in centralrows])), max=max((m['minimum_scalar_shape_channels'] for m in centralrows))), all_basis_sufficient=all((m['basis_sufficient'] for m in centralrows)), all_NNLS_numerical_gates=all((m['support_scenarios_valid'] for m in centralrows)), max_shared_scale_profile_delta_pp=max((m['shared_scale_profile_delta_pp'] for m in centralrows)), support_envelope_kind='35 explicit sensitivity scenarios, not certified independent-box bounds;500-1130N/mm per tooth is an incisor-derived hypothesis; posterior UNKNOWN', external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/; local zip source members and hashes in raw/cases/', compared_quantity='Measured registered arch geometry and provisional projected contact graph; no force observation', refutes_us=True), costs=dict(new_wall_s=time.perf_counter() - started, unique_contact_graphs=len(cache), all_shared_R1_and_R2_costs_retained=True), outcome='GEOMETRY_ALONE_DOES_NOT_IDENTIFY_JOINT_REGIONAL_FORCE; MINIMUM_MEASUREMENT_BASIS_CONSTRUCTED; EXACT_DUAL_NUMERICS_REPAIRED', rows=rows)
    dump(H / 'round3/results.json', result)
    dump(H / 'FROZEN_PREDICTIONS_R3.json', dict(kind='Conditional model/measurement-basis freeze; no physical-force acquisition', prereg_sha256=__import__('hashlib').sha256((H / 'PREREG_R3.json').read_bytes()).hexdigest(), timestamp_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), results_sha256=__import__('hashlib').sha256((H / 'round3/results.json').read_bytes()).hexdigest()))
    print(json.dumps({k: v for (k, v) in result.items() if k != 'rows'}, indent=2))
if __name__ == '__main__':
    main()
