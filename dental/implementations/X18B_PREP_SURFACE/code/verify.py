"""Independent artifact and live-source checks; leaves first results frozen."""
from common import *
from scipy.linalg import qr
from pressure import peak
from surface import read_stl, crop, minimum
from scipy.spatial import cKDTree

def run():
    checks = []
    start = time.perf_counter()
    for tag in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']:
        verify_freeze(tag)
        checks.append(dict(check='Frozen hash chain ' + tag, passed=True))
    result = json.loads((H / 'rounds/R1.json').read_text())
    cur = None
    checked_min = set()
    for r in result['rows']:
        if r['status'] != 'SCORED':
            continue
        if cur != r['case']:
            (data, _) = pair(r['case'])
            cur = r['case']
        orig = data[r['jaw']]['tri'][data[r['jaw']]['owner'] == r['fdi']]
        ex = read_stl(r['external_stl'])
        err = float(np.max(np.abs(orig - ex)))
        checks.append(dict(check=f"Live exterior {r['case']}/{r['fdi']}", passed=err <= 1e-06, max_error_mm=err, injected_shift_rejected=float(np.max(np.abs(orig - (ex + np.array([0, 0, 0.2]))))) > 1e-06))
        z = np.load(r['array_file']['path'])
        (h, _) = query_height(orig, z['xy'], r['jaw'] == 'upper')
        expected = z['source_z']
        finite = np.isfinite(h)
        match = np.array_equal(finite, np.isfinite(expected)) and np.max(np.abs(h[finite] - expected[finite])) <= 1e-09
        checks.append(dict(check=f"Independent original envelope {r['case']}/{r['fdi']}", passed=bool(match)))
        key = (r['case'], r['jaw'])
        if key not in checked_min:
            opp = data['upper' if r['jaw'] == 'lower' else 'lower']['tri']
            opp = crop(opp, z['xy'].min(0), z['xy'].max(0))
            (U, L) = (orig, opp) if r['jaw'] == 'upper' else (opp, orig)
            (fast, full, _) = minimum(U, L)
            checks.append(dict(check=f"Live full-pair minimum {r['case']}/{r['fdi']}", passed=abs(full['minimum_gap_mm'] - r['continuous']['minimum_gap_mm']) <= 1e-07 and abs(fast['minimum_gap_mm'] - full['minimum_gap_mm']) <= 1e-07))
            checked_min.add(key)
    r2 = json.loads((H / 'rounds/R2.json').read_text())
    r3 = json.loads((H / 'rounds/R3.json').read_text())
    for r in r3['rows']:
        op = np.load(r['operator_file'])
        M = op['measurement_operator']
        obs = op['hypothetical_observation']
        (_, _, pivot) = qr(M.T, pivoting=True, mode='economic')
        n = M.shape[1]
        direct = np.linalg.solve(M[pivot[:n]], obs[pivot[:n]])
        prior = next((s for s in r2['stress_results'] if s['case'] == r['case'] and s['fdi'] == r['fdi']))
        basis = np.load(prior['basis']['path'])
        area = basis['area_shares']
        err = float(np.max(np.abs(direct - area)))
        checks.append(dict(check=f"Independent complete channel solve {r['case']}/{r['fdi']}", passed=err <= 1e-09, error=err))
        N = op['nullspace']
        if N.shape[1]:
            direction = N[:, 0]
            allowed = min([area[i] / abs(direction[i]) for i in range(n) if abs(direction[i]) > 1e-12]) * 0.25
            a1 = area + allowed * direction
            a2 = area - allowed * direction
            momenterr = float(np.max(np.abs(op['G'] @ (a1 - a2))))
            tensorerr = float(np.max(np.abs(np.einsum('j,jea->ea', a1 - a2, basis['stress_tensors_MPa']))))
            checks.append(dict(check=f"Same wrench, different regional stress {r['case']}/{r['fdi']}", passed=momenterr <= 1e-09 and tensorerr > 1e-08, moment_difference=momenterr, stress_tensor_difference_MPa=tensorerr))
    r5 = json.loads((H / 'rounds/R5.json').read_text())
    tri = read_stl(r5['preparation_stl']['path'])
    (v, inv) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
    f = inv.reshape(-1, 3)
    e = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    ct = np.unique(e, axis=0, return_counts=True)[1]
    checks.append(dict(check='Exported preparation STL edge incidence', passed=bool(np.all(ct == 2)), nonmanifold_edges=int(np.sum(ct != 2))))
    allok = all((c['passed'] for c in checks))
    dump(H / 'raw/VERIFY.json', dict(status='PASS' if allok else 'FAIL', checks=checks, check_count=len(checks), wall_s=time.perf_counter() - start, original_results_preserved=True))
    print(json.dumps(dict(status='PASS' if allok else 'FAIL', checks=len(checks), wall_s=time.perf_counter() - start), indent=2))
    assert allok
if __name__ == '__main__':
    run()
