from common import *
from pressure import peak
from itertools import combinations
from scipy.linalg import null_space, qr

def run():
    pr = json.loads((H / 'PREREG_R3.json').read_text())
    r2 = json.loads((H / 'rounds/R2.json').read_text())
    original = json.loads((X18 / 'raw/RESULTS_R1_ROWS.json').read_text())
    rows = []
    files = []
    start = time.perf_counter()
    for r in r2['stress_results']:
        z = np.load(r['basis']['path'])
        basis = z['stress_tensors_MPa']
        area = z['area_shares']
        n = len(area)
        rr = next((s for s in original if s.get('status') == 'SCORED' and s['case'] == r['case'] and (s['fdi'] == r['fdi'])))
        patches = rr['arms']['0.1']['original']['patches']
        xy = np.array([p['centroid_xyz_mm'][:2] for p in patches])
        xy -= xy.mean(0)
        G = np.r_[np.ones((1, n)), xy.T]
        rank = int(np.linalg.matrix_rank(G))
        N = null_space(G)
        w = G @ area
        (_, _, rows_pivot) = qr(G.T, pivoting=True, mode='economic')
        E = G[rows_pivot[:rank]]
        b = w[rows_pivot[:rank]]
        verts = []
        for cols in combinations(range(n), rank):
            A = E[:, cols]
            if np.linalg.matrix_rank(A) < rank:
                continue
            xx = np.linalg.solve(A, b)
            if xx.min() >= -1e-10:
                alpha = np.zeros(n)
                alpha[list(cols)] = np.maximum(xx, 0)
                alpha /= alpha.sum()
                if np.max(np.abs(G @ alpha - w)) <= pr['metrics']['moment_parity_max']:
                    verts.append(alpha)
        assert verts
        verts = np.array(verts)
        vertex_peaks = np.array([peak(np.einsum('j,jea->ea', v, basis)) for v in verts])
        conditional_upper = float(vertex_peaks.max())
        area_peak = peak(np.einsum('j,jea->ea', area, basis))
        old_upper = r['old_sharp_upper_MPa']
        (_, _, pivot) = qr(N.T, pivoting=True, mode='economic') if N.shape[1] else (None, None, np.array([], int))
        selected = pivot[:N.shape[1]]
        M = np.r_[G, np.eye(n)[selected]]
        obs = M @ area
        recovered = np.linalg.lstsq(M, obs, rcond=None)[0]
        err = float(np.max(np.abs(recovered - area)))
        badobs = obs.copy()
        badobs[0] *= 2
        badrecovered = np.linalg.lstsq(M, badobs, rcond=None)[0]
        badres = float(np.max(np.abs(M @ badrecovered - badobs)))
        badagreement = float(np.max(np.abs(badrecovered - area)))
        if badres <= pr['metrics']['channel_roundtrip_max']:
            fault_rejected = badagreement > pr['metrics']['channel_roundtrip_max']
            fault_type = 'Reference roundtrip rejects injected resultant; independent redundant loadcell needed for real fault detection'
        else:
            fault_rejected = True
            fault_type = 'Redundant wrench constraint detects injected resultant'
        path = D / f"{r['case']}_{r['fdi']}_calibration_operator.npz"
        np.savez_compressed(path, G=G, nullspace=N, selected_patch_integrals=selected, measurement_operator=M, hypothetical_observation=obs, feasible_vertices=verts, vertex_peaks_MPa=vertex_peaks)
        files.append(artifact(path))
        rows.append(dict(case=r['case'], fdi=r['fdi'], resolution='PER_SURFACE_REGION', patches=n, wrench_rank=rank, unresolved_shape_dimensions=n - rank, extra_patch_integrals_with_known_wrench=n - rank, patch_integrals_with_force_only=n - 1, selected_patch_integrals=selected.tolist(), conditional_moment_upper_MPa_at100N=conditional_upper, unconstrained_simplex_upper_MPa_at100N=old_upper, area_weighted_MPa_at100N=area_peak, conditional_shape_ratio=conditional_upper / area_peak, conditional_upper_reduction_factor=old_upper / conditional_upper, complete_channel_recovery_error=err, injected_resultant_reference_rejection=fault_rejected, injected_resultant_residual=badres, fault_scope=fault_type, physical_measurement_status='NOT_RUN; CoP and patch fractions are hypothetical original-area scenario', operator_file=str(path)))
    prediction = H / 'raw/PREDICTIONS_R3.json'
    dump(prediction, rows)
    files.append(artifact(prediction))
    dump(H / 'FROZEN_PREDICTIONS_R3.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R3.json'), files=files, code_sha256={str(p): sha(p) for p in [H / 'code/calibration.py', H / 'code/common.py', H / 'code/pressure.py']}, physical_measurement_status='NOT_RUN'))
    checks = dict(recovery=all((r['complete_channel_recovery_error'] <= pr['metrics']['channel_roundtrip_max'] for r in rows)), injection=all((r['injected_resultant_reference_rejection'] for r in rows)), bound=all((r['conditional_moment_upper_MPa_at100N'] <= r['unconstrained_simplex_upper_MPa_at100N'] * (1 + 1e-09) for r in rows)))
    out = dict(round='R3', claim_type='capability', external_referent=pr['external_referent'], decision='CALIBRATION_OPERATOR_PASS' if all(checks.values()) else 'CALIBRATION_OPERATOR_FAIL', gates=checks, rows=rows, physical_calibration='NOT_RUN', cost=dict(wall_s=time.perf_counter() - start), scope='Conditional acquisition design: published tooth forces supplied neither CoP nor individual patch integrals. No measured stress narrowing claimed.')
    dump(H / 'rounds/R3.json', out)
    (H / 'HANDOFF_R3.md').write_text('R3 finite contact calibration operator saved. Force+two moments leaves n-rank(G) unresolved dimensions; complementary patch integrals identify the fixed model. Center-of-pressure comparisons are hypothetical, not observed. Next construction R4 changes folded normal-offset preparation to implicit Euclidean erosion of an insertion envelope; actual cervical margin and physical contact calibration still missing.\n')
    state('R3_DECIDED', out['decision'], 'Implicit preparation erosion, then package lab demo')
    print(json.dumps(dict(decision=out['decision'], rows=rows), indent=2))
if __name__ == '__main__':
    run()
