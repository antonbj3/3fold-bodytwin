from diagnose import *
from fractions import Fraction as Q

def check_apex(c):
    ts = c['triangles']
    weights = list(map(Q, c['weights']))
    normals = []
    for t in ts:
        (a, b, z) = [[Q(float(v)) for v in p] for p in t]
        u = [b[k] - a[k] for k in range(3)]
        v = [z[k] - a[k] for k in range(3)]
        n = [-(u[1] * v[2] - u[2] * v[1]), -(u[2] * v[0] - u[0] * v[2]), -(u[0] * v[1] - u[1] * v[0])]
        normals.append(n)
    return len(weights) == len(normals) and min(weights) >= 0 and (max(weights) > 0) and all((sum((y * n[k] for (y, n) in zip(weights, normals))) == 0 for k in range(3))) and all((n == list(map(Q, s)) for (n, s) in zip(normals, c['raw_normals'])))

def run():
    st = time.perf_counter()
    x = read(R / 'RESULTS_A2.json')
    controls = []
    for row in x['rows']:
        c = row['apex']['certificate']
        valid = check_apex(c)
        bad = json.loads(json.dumps(c))
        bad['weights'][0] = str(Q(bad['weights'][0]) + 1)
        badn = json.loads(json.dumps(c))
        badn['raw_normals'][0][0] = str(Q(badn['raw_normals'][0][0]) + 1)
        w = row['native_wall_witness']['exact']
        (a, b) = np.array(w['triangles'])
        new = exact_pair(a, b)
        goodwall = new['distance_squared_mm2'] == w['distance_squared_mm2'] == '0'
        lifted = b + np.array([0.0, 0.0, 100.0])
        shifted = exact_pair(a, lifted)
        controls.append(dict(key=row['key'], exact_positive_dependence=valid, bad_multiplier_rejected=not check_apex(bad), forged_normal_rejected=not check_apex(badn), original_wall_intersection_replayed=goodwall, separated_triangle_not_reported_intersecting=not shifted['fails_0p5mm']))
    tri = np.array([[[-2.0, -2.0, 0], [2.0, -2.0, 0], [2.0, 2.0, 0]], [[-2.0, -2.0, 0], [2.0, 2.0, 0], [-2.0, 2.0, 0]]])
    p = np.array([[0.0, 0.0, 0.05], [0.0, 0.0, 50.0], [0.0, 0.0, -0.05]])
    (q, d, _) = fast_nearest(tri, p)
    signed = p[:, 2]
    film_gate = lambda dd, ss: bool(abs(dd - 0.05) <= 0.01 and ss > 0)
    films = dict(correct=film_gate(d[0], signed[0]), wrong_units_rejected=not film_gate(d[1], signed[1]), wrong_side_rejected=not film_gate(d[2], signed[2]), actual_distances_mm=d.tolist())
    box = trimesh.creation.box()
    broken = box.copy()
    broken.update_faces(np.arange(len(box.faces) - 1))
    topology = dict(closed_box=bool(box.is_watertight), deleted_facet_rejected=not broken.is_watertight)
    a = np.array([0.625, 0.875])
    b = np.array([0.375, 1.125])
    suff = dict(identity_error_mm=float(a.mean() - b.mean()), summary_mean_mm=[float(a.mean()), float(b.mean())], downstream_min_mm=[float(a.min()), float(b.min())], downstream_difference_mm=float(a.min() - b.min()), wall_decision=[bool(a.min() >= 0.5), bool(b.min() >= 0.5)], minimal_extension='Minimum for this fixed decision; retain spatial field for interventions', resolution=dict(summary='PER_TOOTH', local='PER_SURFACE_REGION'))
    allpass = all((all((v for (k, v) in c.items() if k != 'key')) for c in controls)) and all((v for (k, v) in films.items() if k != 'actual_distances_mm')) and all(topology.values()) and (suff['identity_error_mm'] == 0) and (suff['wall_decision'] == [True, False])
    out = dict(all_pass=allpass, exact_controls=controls, film_controls=films, topology_controls=topology, sufficiency=suff, seconds=time.perf_counter() - st, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(R / 'raw/VALIDATION.json', out)
    assert allpass
    print('All exact, corruption, units, sign and topology controls pass', flush=True)
    return out
if __name__ == '__main__':
    run()
