from joint_solid import *
from score_joint import intersections

def run():
    pr = R / 'PREREG_F.json'
    assert not pr.exists()
    dump(pr, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Decide whether the remaining closure failures arise in the artificial patch or immutable measured triangles', obstacle='Protected repair reports success but all13/11/9 self-intersections remain', changed_operation='Classify exact CGAL triangle pairs by native/closure provenance and rerun native patch alone; exact Fraction intersections on nonadjacent pairs', consumer='Decide whether same-source complete valid shell is possible with100% native-facet retention', metrics=dict(max_native_self_intersections=0, exact_intersection_distance_squared_mm2='0', coordinate_replacement_allowed=False), strongest_equally_informed_control='Closed-support versus native-only CGAL; exact triangle-pair witness; clean and deliberately intersecting triangle controls', falsifiers=['All bad pairs involve only artificial closure', 'Native-only mesh has0 intersections', 'Exact rational witness separated'], external_referent=dict(kind='published_code', locator='https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__intersection__grp.html', compared_quantity='Exact-predicate native triangle self-intersection, replayed on local published dataset geometry', refutes_us=True), resolution='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', full_cost=dict(preparation='Already frozen source triangles', fit='None', discovery='3 native tests', validation='Exact witness plus actual corruption control', queries=3, fallback='UNKNOWN if only shared-simplex witness available', historical='UNKNOWN')))
    pr.with_suffix('.sha256').write_text(sha(pr) + '  ' + pr.name + '\n')
    dump(R / 'DECOMPOSITION_F.json', dict(idea='An embedded surface cannot retain two nonadjacent triangles that intersect', equation='Distinct nonincident facets of an embedded2-manifold have disjoint interiors', operation='Exact intersection predicates and source incidence', leaves=[dict(status='DERIVED_UNDER_ASSUMPTIONS', statement='Embedding definition applied to immutable source triangles', stopping_argument='Changing source geometry or accepting intersecting surfaces changes the claim'), dict(status='EXTERNALLY_MEASURED', statement='Native IOS mesh coordinates', stopping_argument='Mesh defect is not anatomical pathology')]))
    idx = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    bad = [r for r in read(R / 'RESULTS_D.json')['rows'] if r['status'] == 'REJECTED']
    rows = []
    for row in bad:
        key = row['key']
        t = load(idx[key]['private_path'])['target']
        (v, inv) = np.unique(t.reshape(-1, 3), axis=0, return_inverse=True)
        f = inv.reshape(-1, 3)
        mm = trimesh.Trimesh(v, f, process=False)
        trimesh.repair.fix_normals(mm, multibody=True)
        (v, f) = (mm.vertices, mm.faces)
        si = intersections(v, f, key + '_F_native')
        pairs = si.get('pairs', [])
        witnesses = []
        for (a, b) in pairs:
            shared = set(f[a]) & set(f[b])
            w = exact_pair(v[f[a]], v[f[b]])
            witnesses.append(dict(face_indices=[a, b], shared_vertex_count=len(shared), exact=w, nonincident_intersection_proved=len(shared) == 0 and w['distance_squared_mm2'] == '0'))
        n_proved = sum((w['nonincident_intersection_proved'] for w in witnesses))
        out = dict(key=key, family=row['family'], native_intersections=si, checked_pairs=witnesses, has_exact_nonincident_obstruction=n_proved > 0, verdict='NO_EMBEDDED_SURFACE_PRESERVING_ALL_NATIVE_FACETS' if n_proved else 'CGAL_INTERSECTION_REQUIRES_INCIDENCE_REVIEW')
        rows.append(out)
        dump(R / 'raw/F_SOURCE_OBSTRUCTIONS.json', rows)
        print(key, si.get('count'), n_proved, out['verdict'], flush=True)
    v = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.25, -1.0], [0.0, 0.25, 1.0], [0.0, -0.5, 0.0]])
    f = np.array([[0, 1, 2], [3, 4, 5]])
    badcheck = intersections(v, f, 'F_injected_crossing')
    v2 = v.copy()
    v2[3:, 0] += 10
    goodcheck = intersections(v2, f, 'F_clean_separated')
    controls = dict(crossing_rejected=badcheck.get('count', 0) > 0, clean_accepted=goodcheck.get('count') == 0)
    assert all(controls.values())
    dump(R / 'RESULTS_F.json', dict(claim_type='capability', rows=rows, controls=controls, external_referent=read(pr)['external_referent']))
if __name__ == '__main__':
    run()
