from local import *
import geometry as g

def run():
    g.closest = fast_nearest
    scorer.closest = fast_nearest
    (a, b) = first()
    m = npz(D6 / 'D_exterior.npz')
    old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
    p = npz(b['public_path'])
    target = npz(b['private_path'])['target']
    (v, f) = (m['vertices'], m['faces'])
    mm = trimesh.Trimesh(v, f, process=False)
    ll = loops(mm)
    shape = scorer.metrics(v[f], target)
    wall = scorer.wall(v[f], old['vertices'][old['faces'][old['roles'] == 1]])
    dd = cKDTree(v).query(p['margin_curve'], workers=1)[0]
    out = dict(shape=shape, wall_against_R5_inner=wall, boundary_loops=[len(x) for x in ll], original_margin_knots_max_distance_mm=float(dd.max()), original_margin_exact=bool((dd == 0).all()), mesh_winding=mm.is_winding_consistent, components=len(mm.split(only_watertight=False)), scope='D exterior only, no complete crown')
    save(R6 / 'RESULTS_D.json', out)
    print(clean(out))
if __name__ == '__main__':
    run()
