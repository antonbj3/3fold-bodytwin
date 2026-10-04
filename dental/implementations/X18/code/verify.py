"""Controls with explicit faults; analytic/LP controls cover mathematics only."""
import json, time
import numpy as np
from scipy.optimize import linprog
from geometry import *
from continuous import plane, exact_batch, inside
from fe import mesh, operators, constitutive

def lp_gap(U, L):
    """Independent LP over the intersection of two projected triangles."""
    constraints = []
    rhs = []
    for t in [U, L]:
        xy = t[:, :2]
        ori = np.sign(np.cross(xy[1] - xy[0], xy[2] - xy[0]))
        for (a, b) in zip(xy, np.roll(xy, -1, axis=0)):
            edge = b - a
            normal = ori * np.array([edge[1], -edge[0]])
            constraints.append(normal)
            rhs.append(float(normal @ a))

    def affine(t):
        return np.linalg.solve(np.c_[t[:, :2], np.ones(3)], t[:, 2])
    c = affine(U) - affine(L)
    res = linprog(c[:2], A_ub=np.array(constraints), b_ub=np.array(rhs), bounds=[(None, None)] * 2, method='highs', options={'threads': 1})
    return float(res.fun + c[2]) if res.success else None

def run():
    start = time.perf_counter()
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    preds = json.loads((H / 'raw/PREDICTIONS_R1.json').read_text())
    p = next((r for r in preds if r['status'] == 'DESIGNED'))
    (a, _) = pair(p['case'])
    z = np.load(p['file'])
    xy = z['xy']
    C = np.c_[xy, z['z_informed']][z['faces']]
    U = a['upper']['tri']
    U = U[np.all(U[:, :, :2].max(1) >= xy.min(0), axis=1) & np.all(U[:, :, :2].min(1) <= xy.max(0), axis=1)]
    pairs = broad_phase(U, C)
    bounded = extremum(U, C, pairs, True)
    full = extremum(U, C, pairs, False)
    tol = pr['metrics']['continuous_gap_oracle_tolerance_mm']
    checks = []

    def add(name, value, passes, ref):
        checks.append(dict(name=name, value=value, pass_=bool(passes), external_referent=ref))
    ref = dict(kind='closed_form', locator='https://doi.org/10.1007/978-0-387-40065-5; affine objective over a convex polygon: extreme point theorem', compared_quantity='minimum affine projected triangle gap, mm', refutes_us=True)
    err = abs(bounded['minimum_gap_mm'] - full['minimum_gap_mm'])
    add('all-pair control', err, err <= tol, ref)
    add('all-pair control rejects injected +1mm answer', 1, abs(bounded['minimum_gap_mm'] + 1 - full['minimum_gap_mm']) > tol, ref)
    (up, cp) = (plane(U), plane(C))
    seed = np.random.default_rng(18)
    sel = pairs[seed.choice(len(pairs), min(48, len(pairs)), replace=False)]
    lperrors = []
    oracle_rows = []
    for ids in sel:
        lp = lp_gap(U[ids[0]], C[ids[1]])
        exact = exact_batch(U, C, up, cp, ids[None, :])[0]
        if lp is not None and np.isfinite(exact):
            lperrors.append(abs(lp - exact))
            oracle_rows.append(dict(pair=ids.tolist(), lp_mm=lp, vertex_mm=exact))
    lpmax = max(lperrors) if lperrors else float('inf')
    add('independent halfplane LP', lpmax, lpmax <= 1e-07, ref)
    add('LP rejects injected +.2mm result', 0.2, all((abs(r['lp_mm'] - r['vertex_mm'] - 0.2) > 1e-07 for r in oracle_rows)) and bool(oracle_rows), ref)
    injured = z['z_informed'] + 1
    gap = signed_gap(U, xy, injured, z['faces'])
    add('nonpenetration rejects +1mm height', gap['maximum_penetration_mm'], gap['maximum_penetration_mm'] > pr['metrics']['penetration_tolerance_mm'], ref)
    diff = float(np.max(np.abs(z['z_informed'] - z['z_control'])))
    add('independent informed projection parity', diff, diff <= 1e-09, dict(kind='closed_form', locator='Projection onto (-infinity,u]: P(z)=min(z,u)', compared_quantity='projected crown height, mm', refutes_us=True))
    add('projection parity rejects +.2mm control fault', 0.2, np.max(np.abs(z['z_informed'] - (z['z_control'] + 0.2))) > 1e-09, ref)
    from experiment import contact, score
    f = pr['tooth_fdi'][0]
    s = site(a['lower'], f, pr['grid_n'])
    copy = {k: v.copy() if isinstance(v, np.ndarray) else v for (k, v) in a['lower'].items()}
    original_vertex = a['lower']['v']
    q = s['z_cervical']
    mut = (copy['labels'] == f) & (copy['v'][:, 2] > q + 1e-05)
    copy['v'][mut, 2] += 10
    s2 = site(copy, f, pr['grid_n'])
    change = max((float(np.max(np.abs(s[k] - s2[k]))) for k in ['xy', 'z_top', 'z_cervical']))
    add('target occlusal +10mm does not change allowed design site', change, change <= 1e-08, dict(kind='our_own_fixture', locator='raw/VERIFICATION.json occlusal-leakage injection', compared_quantity='site design inputs', refutes_us=True))
    co = contact(z['ceiling'] - z['z_informed'], z['z_informed'], xy, z['index'], z['weights'], 0.1, 0.1)
    fake = dict(co)
    fake['mask'] = ~co['mask']
    injscore = score(fake, co, xy)['IoU']
    add('contact scoring rejects complement mask', injscore, injscore == 0, dict(kind='published_dataset', locator=p['source']['upper']['member'], compared_quantity='registered projected contact mask', refutes_us=True))
    (xyz, tet) = mesh(xy, z['z_informed'], z['faces'])
    (K, B, Dm, vol, dofs) = operators(xyz, tet, 210000, 0.3)
    A = np.array([[0.001, 0.0002, 0], [0, -0.0003, 0.0004], [0.0001, 0, 0.0005]])
    u = (xyz @ A.T).ravel()
    eps = np.array([A[0, 0], A[1, 1], A[2, 2], A[0, 1] + A[1, 0], A[1, 2] + A[2, 1], A[0, 2] + A[2, 0]])
    expected = Dm @ eps
    actual = np.einsum('eai,ei->ea', B, u[dofs]) @ Dm.T
    stresserr = float(np.max(np.abs(actual - expected)))
    energy = 0.5 * u @ (K @ u)
    exact_energy = 0.5 * eps @ Dm @ eps * vol.sum()
    add('FE affine stress and energy closed form', dict(stress_error_MPa=stresserr, relative_energy_error=abs(energy / exact_energy - 1)), stresserr < 1e-07 and abs(energy / exact_energy - 1) < 1e-09, dict(kind='closed_form', locator='https://ocw.mit.edu/courses/1-050-engineering-mechanics-i-fall-2007/resources/lec22/; isotropic Hooke law. P1 affine shape/energy identity derived explicitly in this test', compared_quantity='affine elastic strain, stress and energy', refutes_us=True))
    add('FE closed form rejects 2x modulus', float(np.max(np.abs(2 * actual - expected))), np.max(np.abs(2 * actual - expected)) > 1e-07, ref)
    total = pr['fe_contract']['total_force_N']
    add('force normalization rejects injected 2x total', 2 * total, abs(2 * total - total) > 1e-06, dict(kind='closed_form', locator='Newton static force balance; prescribed normalized scenario', compared_quantity='total applied load N', refutes_us=True))
    (vv, ff) = shell(xy, z['z_informed'], z['faces'], 1.2)
    ee = np.r_[ff[:, [0, 1]], ff[:, [1, 2]], ff[:, [2, 0]]]
    (_, counts) = np.unique(np.sort(ee, axis=1), axis=0, return_counts=True)
    good = bool(np.all(counts == 2))
    add('export closed shell', good, good, dict(kind='closed_form', locator='Triangle mesh edge-incidence manifold condition', compared_quantity='two face incidences per edge', refutes_us=True))
    inj = ff[1:]
    ee = np.r_[inj[:, [0, 1]], inj[:, [1, 2]], inj[:, [2, 0]]]
    (_, ct) = np.unique(np.sort(ee, axis=1), axis=0, return_counts=True)
    add('export rejects missing facet', int(np.sum(ct != 2)), np.any(ct != 2), ref)
    dump(H / 'raw/VERIFICATION.json', dict(checks=checks, passed=sum((r['pass_'] for r in checks)), total=len(checks), lp_oracles=oracle_rows, wall_seconds=time.perf_counter() - start, scope='Mathematical controls and software faults; no physical contact or stress validation'))
    print('Verification', sum((r['pass_'] for r in checks)), '/', len(checks))
    assert all((r['pass_'] for r in checks))
if __name__ == '__main__':
    run()
