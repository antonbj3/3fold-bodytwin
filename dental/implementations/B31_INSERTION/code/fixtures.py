import json, copy, pathlib, numpy as np
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def fixture():
    A = np.array([[0.0, -0.25, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.5]])
    z = [-1.0, 0.0, 1.5, 2.0, 2.5, 5.0]
    safe = []
    bulge = []
    for k in range(len(z) - 1):
        a = np.array([[0.1, -1, z[k]], [0.1, 1, z[k]], [0.1, 1, z[k + 1]], [0.1, -1, z[k + 1]]])
        b = a.copy()
        b[b[:, 2] == 2.0, 0] = -0.1
        for ids in [[0, 1, 2], [0, 2, 3]]:
            safe.append(a[ids])
            bulge.append(b[ids])
    return (A, np.array(safe), np.array(bulge), np.array([0, 0, 4.0]))

def run(outdir=R / 'raw'):
    outdir.mkdir(exist_ok=True, parents=True)
    (A, B, C, D) = fixture()
    import trimesh

    def endpoint_minima(triangles):
        mesh = trimesh.Trimesh(triangles.reshape(-1, 3), np.arange(3 * len(triangles)).reshape(-1, 3), process=False)
        return [float(trimesh.proximity.closest_point(mesh, A + shift)[1].min()) for shift in [D * 0, D]]
    endpoint_distances = [endpoint_minima(B), endpoint_minima(C)]
    summary_safe = min(endpoint_distances[0])
    summary_bulge = min(endpoint_distances[1])
    assert endpoint_distances[0] == endpoint_distances[1]
    safe = check_surface(A[None], Obstacle(B), D)
    bad = check_surface(A[None], Obstacle(C), D)
    final_a = check_surface(A[None], Obstacle(B), D * 0)
    final_b = check_surface(A[None], Obstacle(C), D * 0)
    assert summary_safe == summary_bulge and safe['status'] == 'SURFACE_CLEAR_CERTIFIED' and (bad['status'] == 'COLLISION')
    assert final_a['status'] == final_b['status'] == 'SURFACE_CLEAR_CERTIFIED'
    w = bad['witness']
    faults = {}
    x = copy.deepcopy(w)
    x['weights'][0] = '-1'
    faults['negative_weight'] = not validate_witness(np.array(w['source_triangle_mm']), np.array(w['obstacle_triangle_mm']), D, x)
    x = copy.deepcopy(w)
    x['weights'][6] = '0'
    faults['wrong_pose'] = not validate_witness(np.array(w['source_triangle_mm']), np.array(w['obstacle_triangle_mm']), D, x)
    x = copy.deepcopy(w)
    x['weights'][3] = '17'
    faults['wrong_point'] = not validate_witness(np.array(w['source_triangle_mm']), np.array(w['obstacle_triangle_mm']), D, x)
    x = copy.deepcopy(w)
    x['point_mm'][0] += 1.0
    faults['wrong_display_point'] = not validate_witness(np.array(w['source_triangle_mm']), np.array(w['obstacle_triangle_mm']), D, x)
    x = copy.deepcopy(w)
    x['s_exact'] = '0'
    faults['wrong_display_pose'] = not validate_witness(np.array(w['source_triangle_mm']), np.array(w['obstacle_triangle_mm']), D, x)
    faults['omitted_obstacle_abstains'] = check_surface(A[None], Obstacle(np.empty((0, 3, 3))), D)['status'] == 'UNKNOWN'
    controls = []
    for (name, obs) in [('safe', B), ('bulged', C)]:
        cs = [independent_lp(A, t, D) for t in obs]
        hit = any((x['feasible'] for x in cs))
        controls.append(dict(state=name, any_intersection=hit, pair_controls=cs))
        assert hit == (name == 'bulged')
    assert all(faults.values())
    r = dict(claim_type='capability', kind='our_own_fixture', locator='code/fixtures.py: fixture explicit triangles', compared_quantity='closed-form endpoint gap and whole-translation intersection', resolution='PER_POINT', endpoint_summary_mm=[summary_safe, summary_bulge], measured_endpoint_vertex_minima_mm=endpoint_distances, endpoint_distance_scope='All fixture source vertices at both endpoints; plane lower bound x=0.1 proves the true minimum here, not a general vertex-sampling certificate', identity_error=summary_safe - summary_bulge, bulge_mm=0.2, downstream_path=[safe['status'], bad['status']], minimum_extension='retain spatial triangle region and rigid pose/whole-path intersection; endpoint minimum alone insufficient', safe=safe, bulged=bad, endpoint_checks=[final_a, final_b], independent_controls=controls, faults=faults)
    (outdir / 'SUFFICIENCY_AND_FAULTS.json').write_text(json.dumps(r, indent=2) + '\n')
    print('FIXTURE endpoint identity0, intermediate collision detected; all6faults rejected')
if __name__ == '__main__':
    run()
