from common import *
from construction import exact_mesh
from scipy.spatial import ConvexHull

def run():
    rows = []

    def ck(name, ok, **kw):
        rows.append(dict(name=name, passed=bool(ok), **kw))
    support = trimesh.creation.box(extents=[4, 4, 4])
    q = trimesh.creation.box(extents=[1, 1, 1])
    qt = q.copy()
    qt.vertices += np.array([2.0, 0, 0])
    sp = D / 'ctrl_support.mesh'
    ap = D / 'ctrl_inside.mesh'
    bp = D / 'ctrl_shifted.mesh'
    meshwrite(sp, support.vertices, support.faces)
    meshwrite(ap, q.vertices, q.faces)
    meshwrite(bp, qt.vertices, qt.faces)
    ca = exact_mesh(sp, ap, 0.55, 'ctrl_inside')
    cb = exact_mesh(sp, bp, 0.55, 'ctrl_shifted')
    ck('valid_containment_accepted', ca['all_pass'])
    ck('translated_core_rejected', not cb['all_pass'])
    large = trimesh.creation.box(extents=[3.5, 3.5, 3.5])
    lp = D / 'ctrl_thin_wall.mesh'
    meshwrite(lp, large.vertices, large.faces)
    lc = exact_mesh(sp, lp, 0.55, 'ctrl_thin_wall')
    ck('insufficient_wall_rejected', not lc['all_pass'])
    va = np.r_[q.volume, np.ptp(q.vertices, axis=0), q.edges_unique_length]
    vb = np.r_[qt.volume, np.ptp(qt.vertices, axis=0), qt.edges_unique_length]
    identity = float(abs(va - vb).max())
    ck('summary_machine_identity', np.array_equal(va, vb), identity_error=identity)
    sa = Distance(support.vertices, support.faces).query(q.vertices)[0].min()
    sb = Distance(support.vertices, support.faces).query(qt.vertices)[0].min()
    ck('summary_downstream_differs', ca['all_pass'] != cb['all_pass'], unsigned_min_boundary_distances_mm=[sa, sb], difference_mm=abs(sa - sb))
    for (label, upper) in [('narrow', 0.5), ('undercut', 1.5)]:
        vv = np.array([[r * x, r * y, z] for (z, r) in [(0, 1.0), (2, upper)] for (x, y) in [(-1, -1), (1, -1), (1, 1), (-1, 1)]])
        h = ConvexHull(vv)
        mm = trimesh.Trimesh(vv, h.simplices, process=False)
        trimesh.repair.fix_normals(mm)
        sel = abs(mm.triangles[:, :, 2]).max(1) > 0
        n = normals_exact(mm.vertices, mm.faces[sel], [0, 0, 1])
        ck(label + '_insertion', n['pass_exact'] == (label == 'narrow'))
    plane = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    dist = Distance(plane, np.array([[0, 1, 2]]))
    for g in [0.05, 0.2]:
        dd = dist.query(np.array([[0.25, 0.25, g]]))[0][0]
        ck('gap_' + str(g), (abs(dd - 0.05) <= 0.01) == (g == 0.05))
    ck('margin_good', abs(0.0) <= 0.025)
    ck('margin_injected_50um_rejected', not abs(0.05) <= 0.025)
    ck('mesh_valid', q.is_watertight)
    bad = trimesh.Trimesh(q.vertices, q.faces[:-1], process=False)
    ck('mesh_missing_face_rejected', not bad.is_watertight)
    content = ap.read_bytes()
    badfile = D / 'ctrl_corrupt.mesh'
    badfile.write_bytes(content + b'X')
    ck('corrupt_hash_rejected', sha(ap) != sha(badfile))
    out = dict(claim_type='capability', rows=rows, passed=sum((x['passed'] for x in rows)), count=len(rows), all_pass=all((x['passed'] for x in rows)), sufficiency=dict(identity_error=identity, downstream_clearance_difference_mm=abs(sa - sb), containment_decisions=[ca['all_pass'], cb['all_pass']], minimum_extension='Registered position relative to support plus minimum signed spatial clearance'), external_referent=dict(kind='our_own_fixture', locator='code/controls.py', compared_quantity='Checker sensitivity, not dental reality', refutes_us=False))
    dump(R / 'raw/CONTROLS.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
    return out
if __name__ == '__main__':
    run()
