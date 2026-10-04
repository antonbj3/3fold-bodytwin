from dental_release.paths import expand as _release_expand
from common import *
import time, math, itertools
import numpy as np
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
from field_engine.experimental.normal_cone import verify_direction

def peak(t):
    S = np.zeros(t.shape[:-1] + (3, 3))
    S[..., 0, 0] = t[..., 0]
    S[..., 1, 1] = t[..., 1]
    S[..., 2, 2] = t[..., 2]
    S[..., 0, 1] = S[..., 1, 0] = t[..., 3]
    S[..., 1, 2] = S[..., 2, 1] = t[..., 4]
    S[..., 0, 2] = S[..., 2, 0] = t[..., 5]
    return float(max(0, np.linalg.eigvalsh(S)[..., -1].max()))

def polygon(mu):
    circle = [(Q(1), Q(0)), (Q(3, 5), Q(4, 5)), (Q(0), Q(1)), (-Q(3, 5), Q(4, 5)), (-Q(1), Q(0)), (-Q(3, 5), -Q(4, 5)), (Q(0), -Q(1)), (Q(3, 5), -Q(4, 5))]
    inside = [(mu * x, mu * y) for (x, y) in circle]
    outside = []
    for ((a, b), (c, d)) in zip(circle, circle[1:] + circle[:1]):
        det = a * d - b * c
        outside.append((mu * (d - b) / det, mu * (a - c) / det))
    assert all((x * x + y * y == mu * mu for (x, y) in inside))
    assert all((all((a * x + b * y <= mu for (a, b) in circle)) for (x, y) in outside))
    return (circle, inside, outside)

def force_check(normal, tangent, mu, vertical_N, cap_N):
    return normal >= 0 and tangent[0] ** 2 + tangent[1] ** 2 <= mu ** 2 * normal ** 2 and (0 <= vertical_N <= cap_N)

def main():
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_X18_V4.json')
    oldpr = read(BASE / 'LANE_X18_CROWN_ANTAGONIST/PREREG_R3.json')
    contract = oldpr['fe_contract']
    old = read(BASE / 'LANE_X18_CROWN_ANTAGONIST/rounds/R3.json')
    sys.path.insert(0, str(BASE / 'LANE_X18_CROWN_ANTAGONIST/code'))
    fe = module('x72_fe', BASE / 'LANE_X18_CROWN_ANTAGONIST/code/fe.py')
    predictions = read(BASE / 'LANE_X18_CROWN_ANTAGONIST/raw/PREDICTIONS_R1.json')
    originals = read(BASE / 'LANE_X18_CROWN_ANTAGONIST/raw/RESULTS_R1_ROWS.json')
    mu = Q(pr['force_budget_closure']['friction_mu'])
    cap = Q(pr['force_budget_closure']['per_body_vertical_N'])
    (circle, inner, outer) = polygon(mu)
    results = []
    arrays = []
    for s in oldpr['selected_sites']:
        state('X18_RUNNING', 'Per-body vertical20N cap and Coulombmu0.2 frozen', 'FE vector basis ' + str(s))
        source = next((r for r in originals if r.get('status') == 'SCORED' and r['case'] == s['case'] and (r['fdi'] == s['fdi'])))
        pp = next((r for r in predictions if r['case'] == s['case'] and r['fdi'] == s['fdi']))
        zp = record(pp['file'])
        z = np.load(zp)
        rp = record(Path(_release_expand('@DENTAL_WORK_ROOT@/X18_crown_antagonist')) / (str(s['case']) + '_' + str(s['fdi']) + '_reference.npz'))
        ref = np.load(rp)
        zz = ref['original_z'].copy()
        good = np.isfinite(zz)
        if not good.all():
            zz[~good] = zz[good][cKDTree(z['xy'][good]).query(z['xy'][~good])[1]]
        patches = source['arms']['0.1']['original']['patches']
        assert patches
        (xyz, tet) = fe.mesh(z['xy'], zz, z['faces'], contract['thickness_mm'])
        (K, B, Dm, vol, dofs) = fe.operators(xyz, tet, contract['E_MPa'], contract['nu'])
        n = len(z['xy'])
        support = np.flatnonzero(np.sum(z['uv'] ** 2, axis=1) >= 0.7 ** 2)
        fixed = (3 * support[:, None] + np.arange(3)).ravel()
        free = np.setdiff1d(np.arange(len(xyz) * 3), fixed)
        fac = splu(K[free][:, free].tocsc())
        top = 3 * (len(xyz) - n + np.arange(n))
        weights = []
        frames = []
        direction_checks = []
        roof = np.c_[z['xy'], zz]
        tri = roof[z['faces']]
        cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        cross *= np.where(cross[:, 2] < 0, -1, 1)[:, None]
        vn = np.zeros_like(roof)
        for j in range(3):
            np.add.at(vn, z['faces'][:, j], cross)
        vn /= np.linalg.norm(vn, axis=1)[:, None]
        fit_ranks = []
        for patch in patches:
            ix = np.array(patch['point_indices'])
            w = np.zeros(n)
            w[ix] = z['weights'][ix] / z['weights'][ix].sum()
            weights.append(w)
            X = np.c_[z['xy'][ix], np.ones(len(ix))]
            fit_ranks.append(int(np.linalg.matrix_rank(X)))
            normal = np.sum(w[ix, None] * vn[ix], axis=0)
            normal /= np.linalg.norm(normal)
            t1 = np.array([1.0, 0, 0])
            t1 -= normal * (t1 @ normal)
            t1 /= np.linalg.norm(t1)
            t2 = np.cross(normal, t1)
            frames.append(np.array([normal, t1, t2]))
        rational_normals = [[-Q(float(v)) for v in f[0]] for f in frames]
        escape = verify_direction(rational_normals, [0, 0, -1])
        bad_escape = verify_direction(rational_normals, [0, 0, 1])
        assert escape.status == 'YES_STRICT' and bad_escape.status == 'UNKNOWN'

        def tensor(force_field):
            f = np.zeros(len(xyz) * 3)
            for a in range(3):
                f[top + a] = force_field[:, a]
            u = np.zeros_like(f)
            u[free] = fac.solve(f[free])
            sig = np.einsum('eai,ei->ea', B, u[dofs]) @ Dm.T
            balance = np.sum((K @ u - f)[fixed].reshape(-1, 3), axis=0) + force_field.sum(0)
            return (sig, float(np.linalg.norm(balance)), float(np.linalg.norm((K @ u - f)[free])))
        basis = []
        for w in weights:
            bb = []
            for axis in range(3):
                ff = np.zeros((n, 3))
                ff[:, axis] = w
                (ss, _, _) = tensor(ff)
                bb.append(ss)
            basis.append(bb)
        basis = np.array(basis)
        oldbasis = np.load(record(next((r for r in old['rows'] if r['case'] == s['case'] and r['fdi'] == s['fdi']))['basis_file']))['stress_tensors_MPa']
        vertical = -100 * basis[:, 2]
        oldres = float(np.max(np.abs(vertical - oldbasis)) / max(1, np.max(np.abs(oldbasis))))
        assert oldres <= 1e-08
        local_vertical = max((peak(b * float(cap / 100)) for b in vertical))
        directions = {}
        values = {}
        force_vertices = {}
        for (label, poly) in [('tilted', [(Q(0), Q(0))]), ('friction_inner', inner), ('friction_outer', outer)]:
            tensors = []
            forces = []
            ids = []
            for (j, frame) in enumerate(frames):
                for uv in poly:
                    ray = frame[0] + float(uv[0]) * frame[1] + float(uv[1]) * frame[2]
                    assert ray[2] > 0, 'VERTICAL_BUDGET_DOES_NOT_BOUND_THIS_CONE'
                    force = -float(cap) * ray / ray[2]
                    stress = np.einsum('a,aei->ei', force, basis[j])
                    tensors.append(stress)
                    forces.append(force)
                    ids.append(j)
            pv = [peak(v) for v in tensors]
            imax = int(np.argmax(pv))
            values[label] = max(pv)
            directions[label] = dict(patch=ids[imax], force_N=forces[imax], tangent_poly=poly, prediction_MPa=pv[imax])
            force_vertices[label] = dict(forces_N=forces, patch_ids=ids, peak_MPa=pv)
        assert values['friction_inner'] <= values['friction_outer'] * (1 + 1e-12)
        area = np.array([p['area_mm2'] for p in patches])
        alpha = area / area.sum()
        area_vertical = peak(np.einsum('j,jei->ei', alpha, vertical) * float(cap / 100))
        area_tilted = peak(sum((np.einsum('a,aei->ei', -float(cap) * frame[0] / frame[0, 2], basis[j]) * alpha[j] for (j, frame) in enumerate(frames))))
        mixalpha = np.random.default_rng(s['case'] * 100 + s['fdi']).dirichlet(np.ones(len(patches)))
        forces = [-float(cap) * (fr[0] + float(mu) * fr[1]) / (fr[0, 2] + float(mu) * fr[1, 2]) for fr in frames]
        ff = sum((w[:, None] * f[None, :] * v for (w, f, v) in zip(weights, forces, mixalpha)))
        (direct, bal, res) = tensor(ff)
        combo = sum((np.einsum('a,aei->ei', f, basis[j]) * v for (j, (f, v)) in enumerate(zip(forces, mixalpha))))
        scale = max(1, np.max(np.abs(direct)))
        parity = float(np.max(np.abs(direct - combo)) / scale)
        bad_tensor = float(np.max(np.abs(direct - 2 * combo)) / scale)
        assert parity <= 1e-08 and bad_tensor > 1e-08 and (bal <= 1e-06)
        original = next((r for r in old['rows'] if r['case'] == s['case'] and r['fdi'] == s['fdi']))
        ratio = values['friction_outer'] / max(area_tilted, 1e-12)
        decision = 'UNKNOWN_FORCE_SENSITIVE' if ratio > pr['unchanged_decision_ratio'] else 'CONDITIONAL_LOW_SENSITIVITY'
        path = DATA / 'facet_normals' / (str(s['case']) + '_' + str(s['fdi']) + '_vector_basis.npz')
        path.parent.mkdir(exist_ok=True)
        np.savez_compressed(path, stress_world_vector_unit_N=basis, weights=np.array(weights), frames=np.array(frames), tet=tet, xy=z['xy'])
        arrays.append(dict(path=path, sha256=sha(path), bytes=path.stat().st_size))
        v0 = peak(-float(cap) * basis[0, 2])
        v1 = peak(-float(cap) * basis[1, 2])
        summary = [sum([cap] + [Q(0)] * (len(patches) - 1)), sum([Q(0), cap] + [Q(0)] * (len(patches) - 2))]
        assert summary[0] == summary[1]
        row = dict(**s, patches=len(patches), before_peak_MPa=original['sharp_upper_peak_MPa'], local_vertical_peak_MPa=local_vertical, local_tilted_peak_MPa=values['tilted'], coulomb_peak_interval_MPa=[values['friction_inner'], values['friction_outer']], area_vertical_20N_peak_MPa=area_vertical, area_tilted_20N_peak_MPa=area_tilted, ratio_before=original['upper_over_area_peak'], ratio_after_outer=ratio, decision_before=original['status'], decision_after=decision, changed_decision=decision != original['status'], normal_tilt_deg=[math.degrees(math.acos(f[0, 2])) for f in frames], local_normal_cone=dict(status=escape.status, checked_normals=rational_normals, direction=[0, 0, -1], opposite_direction_rejected=True, scope='local signed gap increases when lower crown moves down with antagonist fixed; no full swept insertion or force-budget inference'), worst_loads=directions, local_budget=dict(crown_vertical_N=str(cap), opposing_aggregate_body_vertical_N=str(cap), scene_vertical_capacity_N=100, global_load_equality='NOT_ASSERTED; these eight sites mix two source cases'), direct_mixed_control=dict(relative_tensor_error=parity, bad2x_tensor_rejected=True, force_balance_N=bal, equilibrium_residual_N=res), stored_vertical100N_relative_error=oldres, sufficiency=dict(summary='local total vertical force20N', values=summary, identity_error_N=0, downstream_peak_difference_MPa=abs(v0 - v1), minimal_extension='patch force allocation plus local normal/tangential force; total body force is insufficient'), rigorous_enclosure='Exact geometric inner disk/outer tangent polygon sandwich and convexity proof; FE/eigenvalue/normal-fit floating rounding enclosure MISSING; not continuum stress validation')
        results.append(row)
        print('X18', s['fdi'], 'old', round(row['before_peak_MPa'], 2), 'local/friction', list(map(lambda x: round(x, 2), row['coulomb_peak_interval_MPa'])), 'decision', decision, flush=True)
        row['plane_fit_rank_failures'] = sum((r < 3 for r in fit_ranks))
        row['normal_construction'] = 'area-weighted incident facets then pressure-weighted patch mean; same mesh roof'
    injections = dict(local_cap_violation_rejected=not force_check(Q(21), (0, 0), mu, Q(21), cap), friction_violation_rejected=not force_check(Q(20), (Q(20), 0), mu, cap, cap), normal_direction_violation_rejected=all((r['local_normal_cone']['opposite_direction_rejected'] for r in results)), tensor2x_rejected=all((r['direct_mixed_control']['bad2x_tensor_rejected'] for r in results)))
    assert all(injections.values())
    prediction = dict(frozen_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_X18_V4.json'), rows=results, arrays=arrays, prospective_measurement_status='NOT_RUN; retrospective same-model comparisons only')
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', prediction)
    stored = read(ROOT / 'FROZEN_PREDICTIONS.json')
    assert stored['rows'] == enc(results) and stored['arrays'] == enc(arrays), 'FROZEN_PREDICTION_DRIFT'
    out = dict(claim_type='capability', resolution='PER_SURFACE_REGION', rows=results, injections=injections, arrays=arrays, changed_conclusion=any((r['changed_decision'] for r in results)), external_referent=dict(kind='published_code', locator='Field7abf765:normal_cone.py and original X18 linear roof FE; convexity: https://web.stanford.edu/~boyd/cvxbook/', compared_quantity='exact local normal-cone direction; convex maximum over linear elastic force polytope', refutes_us=True), independent_stress_measurement='MISSING; source Bits2Bites geometry is not stress/load facit', phenomenological_debts=dict(local20N_cap='replace by same-case calibrated per-body force budget; no empirical20N claim', friction_mu02='measure wet material-pair traction ratio at relevant rate and normal pressure', patch_normal='replace fitted one normal per patch by calibrated source facet normals and orientation uncertainty', elastic_roof='measure full crown/preparation/PDL support and force-displacement; quantify FE refinement'), dropout=dict(selected_sites=8, retained=8, rejected=0, body_scene_aggregation_rejected='7case127+1case122 cannot be one arch'), control='direct FE tensor parity on same roof; vertical100N predecessor retained', cost_wall_s=time.perf_counter() - start)
    write(ROOT / 'raw/X18.json', out)
    state('X18_DONE', 'Budget+cone+mixed-FE fault controls pass; physical stress UNKNOWN', 'Compile before/after report and graph feedback')
if __name__ == '__main__':
    main()
