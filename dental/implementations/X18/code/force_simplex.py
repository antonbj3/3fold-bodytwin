"""A force-uncertainty capability on actual original scan proximity patches."""
import time, datetime
from scipy.sparse.linalg import splu
from geometry import *
from fe import mesh, operators

def peak(tensor):
    S = np.zeros(tensor.shape[:-1] + (3, 3))
    S[..., 0, 0] = tensor[..., 0]
    S[..., 1, 1] = tensor[..., 1]
    S[..., 2, 2] = tensor[..., 2]
    S[..., 0, 1] = S[..., 1, 0] = tensor[..., 3]
    S[..., 1, 2] = S[..., 2, 1] = tensor[..., 4]
    S[..., 0, 2] = S[..., 2, 0] = tensor[..., 5]
    return float(max(0, np.max(np.linalg.eigvalsh(S)[..., -1])))

def run():
    start = time.perf_counter()
    pr = json.loads((H / 'PREREG_R3.json').read_text())
    contract = pr['fe_contract']
    preds = json.loads((H / 'raw/PREDICTIONS_R1.json').read_text())
    rows = json.loads((H / 'raw/RESULTS_R1_ROWS.json').read_text())
    out = []
    files = []
    existing = json.loads((H / 'rounds/R3.json').read_text()) if (H / 'FROZEN_PREDICTIONS_R3.json').exists() else None
    for s in pr['selected_sites']:
        rr = next((r for r in rows if r.get('status') == 'SCORED' and r['case'] == s['case'] and (r['fdi'] == s['fdi'])))
        p = next((r for r in preds if r['case'] == s['case'] and r['fdi'] == s['fdi']))
        z = np.load(p['file'])
        ref = np.load(D / (str(s['case']) + '_' + str(s['fdi']) + '_reference.npz'))
        zz = ref['original_z'].copy()
        good = np.isfinite(zz)
        if not good.all():
            zz[~good] = zz[good][cKDTree(z['xy'][good]).query(z['xy'][~good])[1]]
        patches = rr['arms']['0.1']['original']['patches']
        n = len(z['xy'])
        (xyz, tet) = mesh(z['xy'], zz, z['faces'])
        (K, B, Dm, vol, dofs) = operators(xyz, tet, contract['E_MPa'], contract['nu'])
        support = np.flatnonzero(np.sum(z['uv'] ** 2, axis=1) >= 0.7 ** 2)
        fixed = (3 * support[:, None] + np.arange(3)).ravel()
        free = np.setdiff1d(np.arange(len(xyz) * 3), fixed)
        factor = splu(K[free][:, free].tocsc())
        loads = []
        basis = []
        areas = []

        def stress(load):
            force = np.zeros(len(xyz) * 3)
            force[3 * (len(xyz) - n + np.arange(n)) + 2] = -load
            u = np.zeros_like(force)
            u[free] = factor.solve(force[free])
            return np.einsum('eai,ei->ea', B, u[dofs]) @ Dm.T
        for patch in patches:
            ix = np.asarray(patch['point_indices'])
            load = np.zeros(n)
            load[ix] = 100 * z['weights'][ix] / z['weights'][ix].sum()
            loads.append(load)
            basis.append(stress(load))
            areas.append(patch['area_mm2'])
        if not patches:
            out.append(dict(**s, status='UNKNOWN_NO_CONTACT'))
            continue
        basis = np.array(basis)
        loads = np.array(loads)
        shares = np.array(areas) / sum(areas)
        area_tensor = np.einsum('j,jea->ea', shares, basis)
        vertex_peaks = np.array([peak(t) for t in basis])
        upper = float(vertex_peaks.max())
        area_peak = peak(area_tensor)
        alpha = np.random.default_rng(s['case'] * 100 + s['fdi']).dirichlet(np.ones(len(patches)))
        direct = stress(alpha @ loads)
        combined = np.einsum('j,jea->ea', alpha, basis)
        scale = max(1, float(np.max(np.abs(direct))))
        resid = float(np.max(np.abs(direct - combined)) / scale)
        false_resid = float(np.max(np.abs(direct - 2 * combined)) / scale)
        mixpeak = peak(direct)
        flat = (basis[1:] - basis[0]).reshape(max(0, len(basis) - 1), -1)
        rank = int(np.linalg.matrix_rank(flat, tol=1e-08 * np.linalg.norm(flat, 2))) if len(basis) > 1 else 0
        path = D / (str(s['case']) + '_' + str(s['fdi']) + '_force_basis.npz')
        if existing is None:
            np.savez_compressed(path, stress_tensors_MPa=basis, loads_N=loads, area_shares=shares, original_roof_z=zz, xy=z['xy'], tet=tet)
        else:
            saved = np.load(path)
            assert np.allclose(saved['stress_tensors_MPa'], basis, rtol=1e-09, atol=1e-09)
            assert np.array_equal(saved['loads_N'], loads)
        files.append(dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size))
        ratio = upper / max(area_peak, 1e-12)
        out.append(dict(**s, status='UNKNOWN_FORCE_SENSITIVE' if ratio > pr['metrics']['unobservable_ratio_threshold'] else 'CONDITIONAL_LOW_SENSITIVITY', patch_count=len(patches), area_weighted_peak_MPa=area_peak, sharp_upper_peak_MPa=upper, upper_over_area_peak=ratio, vertex_peaks_MPa=vertex_peaks, unit_total_force_N=100, force_sum_max_error_N=float(np.max(np.abs(loads.sum(1) - 100))), full_tensor_shape_rank=rank, patch_force_shape_channels_sufficient=len(patches) - 1, total_force_shape_channels=0, independent_mixed_solve=dict(alpha=alpha, peak_MPa=mixpeak, superposition_relative_residual=resid, within_sharp_bound=mixpeak <= upper * (1 + 1e-10), injected2x_tensor_relative_residual=false_resid, fault_rejected=false_resid > pr['metrics']['superposition_relative_tolerance']), basis_file=str(path)))
        print('Force simplex', s, 'ratio', ratio, flush=True)
    checks = all((r.get('status') == 'UNKNOWN_NO_CONTACT' or (r['independent_mixed_solve']['superposition_relative_residual'] <= pr['metrics']['superposition_relative_tolerance'] and r['independent_mixed_solve']['fault_rejected'] and r['independent_mixed_solve']['within_sharp_bound'] and (r['force_sum_max_error_N'] <= pr['metrics']['force_conservation_N_max'])) for r in out))
    output = dict(round='R3', claim_type='capability', external_referent=pr['external_referent'], rows=out, mathematical_gate='PASS' if checks else 'FAIL', physical_stress_accuracy='UNKNOWN', cost=dict(wall_seconds=time.perf_counter() - start, peak_rss_kib=__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss), scope='Sharp maximum only within fixed linear roof and uniform-pressure observed proximity patches. No pressure/stress measurement. Sufficient force-channel count is not a proved minimum for a single nonlinear peak.')
    if existing is None:
        dump(H / 'rounds/R3.json', output)
        dump(H / 'FROZEN_PREDICTIONS_R3.json', dict(round='R3', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R3.json'), result_sha256=sha(H / 'rounds/R3.json'), files=files, physical_measurement_status='NOT_RUN', tensor_basis_ready_for_future_registered_patch_pressure=True))
    else:
        fr = json.loads((H / 'FROZEN_PREDICTIONS_R3.json').read_text())
        assert sha(H / 'rounds/R3.json') == fr['result_sha256']
        assert all((sha(f['path']) == f['sha256'] for f in fr['files']))
        dump(H / 'raw/REPLAY_R3_VALIDATION.json', dict(status='PASS' if checks else 'FAIL', cost=output['cost'], frozen_prediction_preserved=True, independent_direct_solves=len(out)))
    assert checks
if __name__ == '__main__':
    run()
