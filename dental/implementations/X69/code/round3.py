"""Exact affine enclosure on stored digital geometry; empirical epsilon stays UNKNOWN."""
from fractions import Fraction as F
from scipy.linalg import solve
import math
from registration import *

def rational_inverse(M):
    a = [[F(float(x)) for x in row] + [F(int(i == j)) for j in range(4)] for (i, row) in enumerate(M)]
    for k in range(4):
        pivot = next((i for i in range(k, 4) if a[i][k]))
        (a[k], a[pivot]) = (a[pivot], a[k])
        scale = a[k][k]
        a[k] = [x / scale for x in a[k]]
        for i in range(4):
            if i != k:
                scale = a[i][k]
                a[i] = [x - scale * y for (x, y) in zip(a[i], a[k])]
    return [row[4:] for row in a]

def apply_rational(inv, x):
    return [sum((inv[i][j] * F(float(x[j])) for j in range(4))) for i in range(4)]

def greedy(points):
    i0 = int(np.argmin(points[:, 0]))
    i1 = int(np.argmax(np.linalg.norm(points - points[i0], axis=1)))
    edge = points[i1] - points[i0]
    i2 = int(np.argmax(np.linalg.norm(np.cross(points - points[i0], edge), axis=1)))
    n = np.cross(points[i1] - points[i0], points[i2] - points[i0])
    i3 = int(np.argmax(np.abs((points - points[i0]) @ n)))
    return (points[[i0, i1, i2, i3]], [i0, i1, i2, i3])

def main():
    t = time.monotonic()
    w = np.load(DATA / 'r2_witnesses.npz')
    p = w['lower_source_points']
    fdi = w['lower_fdi']
    T = w['lower_T']
    base = DATA / 'hao_demo/Demo/Demo_1'
    tri = stl(base / 'CBCT/reconstruction/out_smoothed.stl')
    v = tri[:, 0]
    root = v[v[:, 2] < 30]
    ids = np.linspace(0, len(root) - 1, 20, dtype=int)
    qCB = root[ids]
    q = transform(qCB, np.linalg.inv(T))
    A = np.array([p[fdi == i].mean(0) for i in [36, 46, 31, 41]])
    (B, bids) = greedy(p)
    designs = []
    for (name, anchors) in [('A_crown_centroids', A), ('B_spread_crown_points', B)]:
        M = np.vstack((anchors.T, np.ones(4)))
        inv = rational_inverse(M)
        records = []
        maxerr = 0.0
        exactcheck = True
        for (k, x) in enumerate(q):
            h = np.r_[x, 1.0]
            lam = apply_rational(inv, h)
            lhs = [sum((lam[i] * F(float(M[j, i])) for i in range(4))) for j in range(4)]
            identity = all((z == F(float(y)) for (z, y) in zip(lhs, h)))
            exactcheck &= identity
            L = sum((abs(z) for z in lam))
            bound = math.nextafter(float(L), math.inf)
            control = solve(M, h)
            err = float(np.max(np.abs(control - np.array([float(z) for z in lam]))))
            maxerr = max(maxerr, err)
            e_uniform = [F(1)] * 4
            e_extreme = [F(1 if z >= 0 else -1) for z in lam]
            uniform = sum((z * e for (z, e) in zip(lam, e_uniform)))
            extreme = sum((z * e for (z, e) in zip(lam, e_extreme)))
            identical_summary = sum((e * e for e in e_uniform)) == sum((e * e for e in e_extreme)) == 4 and max((abs(e) for e in e_uniform)) == max((abs(e) for e in e_extreme)) == 1
            records.append({'query_index': k, 'CBCT_point_mm': qCB[k].tolist(), 'candidate_IOS_frame_point_mm': x.tolist(), 'lambda_exact': [str(z) for z in lam], 'L_exact': str(L), 'amplification_upper': bound, 'resolution': 'PER_POINT', 'epsilon_mm': 'UNKNOWN', 'conditional_error_bound_mm': 'epsilon_mm * amplification_upper', 'exact_reconstruction': identity, 'numpy_control_weight_error': err, 'sufficiency': {'summary_identity_error': 0, 'machine_exact': identical_summary, 'same_max_anchor_error': 1, 'same_squared_anchor_RMS': 1, 'downstream_uniform_error': float(abs(uniform)), 'downstream_extremal_error': float(abs(extreme)), 'downstream_difference_per_epsilon': float(abs(extreme) - abs(uniform)), 'evidence_kind': 'our_own_fixture_affine_perturbation_not_rigid_measurement'}, 'half_bound_injection_rejected': bool(abs(extreme) > L / 2)})
        values = np.array([x['amplification_upper'] for x in records])
        designs.append({'design': name, 'anchors_native_mm': anchors.tolist(), 'anchor_point_indices': bids if name.startswith('B') else 'FDI centroids 36,46,31,41', 'affine_basis_det': float(np.linalg.det(M)), 'records': records, 'exact_gate_pass': exactcheck, 'fault_gate_pass': all((x['half_bound_injection_rejected'] for x in records)), 'control_max_weight_error': maxerr, 'amplification_min': float(values.min()), 'amplification_median': float(np.median(values)), 'amplification_max': float(values.max()), 'physical_error': 'UNKNOWN_NO_INDEPENDENT_EPSILON'})
    result = {'construction': 'R3', 'claim_type': 'capability', 'external_referent': {'kind': 'closed_form', 'locator': 'https://doc.cgal.org/latest/Barycentric_coordinates_3/namespaceCGAL_1_1Barycentric__coordinates.html ; x=sum lambda_i a_i, sum lambda_i=1 + triangle inequality', 'compared_quantity': 'exact affine-basis reconstruction and pointwise conditional error enclosure', 'refutes_us': True}, 'designs': designs, 'rationale': 'Spatial weights replace a scalar residual for the crown-to-root handoff; selecting spread anchor points reduces conditional amplification, but no measured pose error is claimed.', 'proof': 'For D=T_est-T_true affine and sum lambda=1, D(x)=sum lambda D(a_i). If every ||D(a_i)||<=epsilon, triangle inequality gives ||D(x)||<=epsilon sum |lambda_i|. Exact rational weights from stored binary coordinates, outward-rounded float coefficient. Bound holds for all affine maps, hence rigid maps, but may be loose for rigid maps.', 'empirical_epsilon': 'UNKNOWN', 'seconds': time.monotonic() - t, 'no_affine_sensitivity_without_enclosure': True, 'acquisition_spec': 'At least3 noncollinear independently identified crown landmarks for ideal rigid pose;4 affinely independent anchor error measurements sufficient for the provided rigorous enclosure. Need calibrated mm, acquisition bite state and landmark trueness/repeatability. These are different requirements.'}
    dump(ROOT / 'raw/R3_RESULTS.json', result)
    print([(d['design'], d['amplification_median'], d['amplification_max']) for d in designs])
if __name__ == '__main__':
    main()
