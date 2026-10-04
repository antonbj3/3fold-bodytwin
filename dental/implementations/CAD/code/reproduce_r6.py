import thread_guard
from cadlib import *
from readonly import parents
import trimesh
P = parents()
old = {r['uid']: r for r in read(ROOT / 'raw/R5_PREDICTIONS.json')}
inputs = {r['key']: r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']}
out = []
for prior in read(ROOT / 'RESULTS_R6.json')['rows']:
    r = old[prior['uid']]
    a = np.load(r['mesh_path'])
    dat = np.load(inputs[r['key']]['public_geometry'])
    tri = dat['donor'].copy()
    tri[:, :, 0] *= -1
    tri = tri[:, ::-1]
    center = tri.reshape(-1, 3).mean(0)
    tri[:, :, :2] -= center[:2]
    m = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(tri.size // 3).reshape(-1, 3), process=True)
    fit = r['donor_fit']
    m.vertices[:, :2] *= fit['scale_xy']
    m.vertices = m.vertices @ np.array(fit['R']).T + fit['translation_mm']
    (m, drop, n) = P['r4'].component(m.triangles)
    m = m.slice_plane([0, 0, a['iv'][:, 2].min()], [0, 0, 1], cap=False)
    m.merge_vertices(digits_vertex=10)
    rec = dict(uid=r['uid'], source_component_area_drop=drop, source_components=n, vertices=len(m.vertices), faces=len(m.faces))
    try:
        loops = P['r4'].loops(m)
        rec['boundary_count'] = len(loops)
        rec['boundary_vertex_counts'] = [len(x) for x in loops]
        reason = 'UNKNOWN_DONOR_TOPOLOGY_OR_COMPONENT_LOSS' if drop > 0.05 or len(loops) != 1 else 'UNEXPECTED_PASS'
    except ValueError as e:
        reason = str(e)
    edges = np.sort(m.edges, axis=1)
    (u, counts) = np.unique(edges, axis=0, return_counts=True)
    bd = u[counts == 1]
    degrees = np.bincount(bd.ravel(), minlength=len(m.vertices))
    rec.update(reason=reason, boundary_branch_vertices=int(np.sum((degrees > 0) & (degrees != 2))), nonmanifold_edges=int(np.sum(counts > 2)), same_failure=reason in prior.get('reason', ''))
    if not rec['same_failure']:
        raise AssertionError('R6 source failure changed')
    out.append(rec)
dump(ROOT / 'raw/R6_TOPOLOGY_REPLAY.json', dict(rows=out, all_same_failure=True, resolution='PER_POINT boundary incidence and PER_TOOTH rejection'))
print('R6_REPLAY', len(out), 'same refusals')
