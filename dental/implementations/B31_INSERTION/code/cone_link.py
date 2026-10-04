"""Necessary local virtual-preparation condition linked to whole neighbour sweep.
R4 cavity boundary used as coincident virtual die boundary. No separate physical prep.
"""
import pathlib, json, datetime, hashlib, time
import numpy as np
import warnings
from sweep import Q, rational, linprog
from fixtures import fixture
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze():
    p = R / 'PREREG_R5_CONE_LINK.json'
    p.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Test whether neighbouring-surface-cleared straight axes satisfy the existing virtual intaglio local withdrawal condition', changed_operation='Link full-path axes back to every original R4intaglio facet; exact dyadic normal dot D, independent fixed-D LP', obstacle='An axis that avoids neighbours may immediately close a coincident virtual prep cavity', assumptions='Original winding defines cavity-inward normals; coincident virtual prep boundary is a constitutive scenario, not a measured preparation. Nonnegative local condition necessary only, never sufficient for a nonconvex entire die.', gate='Every original18baseline and2selected axes retained; any negative exact normaldot rejects local condition; LPbooleans must agree; no full virtual/physicalPASS; identical normal summary counterexample required', falsifier='Flip reported negative dot to positive -> exact input-bound verifier rejects', full_cost='20smallfacetsets +20LPs, zero fit; prior whole sweeps reused with hashes; physicalprep/registration UNKNOWN', resolution='PER_POINT', timescale='SIMULTANEOUS', leaves=[dict(status='DERIVED_UNDER_ASSUMPTIONS', equation='n_prep=-cross(b-a,c-a); local necessary n_prep.D>=0', stop='Facet orientation/closure conditional; negative dot for coincident boundary implies local closing for sufficiently small displacement'), dict(status='CONSTITUTIVE_CLOSURE', name='R4 intaglio as coincident virtual prep boundary', stop='No measured distinct prep; no empirical spacer introduced'), dict(status='UNKNOWN', name='complete prep solid and selfintersection-free shell', stop='Local normal inequalities cannot certify entire trajectory')]), indent=2) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def exact_normal_dot(t, D):
    (a, b, c) = [[rational(x) for x in row] for row in t]
    u = [b[j] - a[j] for j in range(3)]
    v = [c[j] - a[j] for j in range(3)]
    n = [-(u[1] * v[2] - u[2] * v[1]), -(u[2] * v[0] - u[0] * v[2]), -(u[0] * v[1] - u[1] * v[0])]
    return sum((n[j] * rational(D[j]) for j in range(3)))

def run(root=R):
    root = pathlib.Path(root)
    tick = time.perf_counter()
    co = json.load(open(root / 'FROZEN_COHORT.json'))
    axes = json.load(open(root / 'rounds/R4/AXES.json'))
    selected = {r['key']: r['selected_translation_mm'] for r in axes['rows'] if r['selected_translation_mm']}
    rows = []
    mismatches = 0
    for rec in co['rows']:
        m = np.load(rec['mesh_path'])
        ids = np.flatnonzero(m['roles'] == 1)
        tri = m['vertices'][m['faces'][ids]]
        directions = [('baseline', [0, 0, 40])] + ([('neighbour_clear', selected[rec['key']])] if rec['key'] in selected else [])
        for (name, D) in directions:
            dots = [exact_normal_dot(t, D) for t in tri]
            k = min(range(len(dots)), key=dots.__getitem__)
            admissible = min(dots) >= 0
            n = -np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                s = linprog([0, 0, 0], A_ub=-n, b_ub=np.zeros(len(n)), bounds=[(x, x) for x in D], method='highs', options={'threads': 4})
            mismatches += int(bool(s.success) != admissible)
            rows.append(dict(key=rec['key'], axis_kind=name, translation_mm=D, local_necessary_condition='SATISFIED_ONLY' if admissible else 'REFUTED_CONSTITUTIVE_VIRTUAL_BOUNDARY', worst_original_face=int(ids[k]), unnormalized_normal_dot_D_exact_mm3=str(dots[k]), independent_LP_feasible=bool(s.success), full_path='UNKNOWN unless already source-neighbour COLLISION; local condition is not whole die clearance', resolution='PER_POINT', timescale='SIMULTANEOUS', reported_sign_flip_rejected=exact_normal_dot(tri[k], D) != -dots[k] if not admissible else None))
    assert mismatches == 0
    (A, B, C, D) = fixture()
    summary = exact_normal_dot(A, D)
    assert summary == 0
    suff = dict(summary_exact=[str(summary), str(summary)], identity_error_exact='0', same_crown_normal_source=True, downstream_path=['SURFACE_CLEAR_CERTIFIED', 'COLLISION'], measured_by='raw/SUFFICIENCY_AND_FAULTS.json', minimum_extension='Full spatial obstacle scene and swept crown, not additional scalar normal margin', kind='our_own_fixture')
    out = dict(claim_type='capability', rows=rows, reference_mismatches=mismatches, sufficiency=suff, seconds=time.perf_counter() - tick, physical='UNKNOWN', cohort_sha256=sha(root / 'FROZEN_COHORT.json'), axis_results_sha256=sha(root / 'rounds/R4/AXES.json'), external_referent={'kind': 'published_code', 'locator': 'https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html', 'compared_quantity': 'fixed-vector facet-normal linear feasibility compared with exact rational dot signs', 'refutes_us': True})
    (root / 'raw/CONE_LINK.json').write_text(json.dumps(out, indent=2) + '\n')
    print('R5CONE_LINK', [(r['key'], r['local_necessary_condition']) for r in rows if r['axis_kind'] == 'neighbour_clear'], flush=True)
    return out
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    a = ap.parse_args()
    freeze() if a.freeze else run()
