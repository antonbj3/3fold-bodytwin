import resource
resource.setrlimit(resource.RLIMIT_AS, (3584 * 1024 ** 2, 3584 * 1024 ** 2))
from common_r3 import *
sys.path.insert(0, str(OLD / 'code'))
import score as inherited
from memory_height import height
inherited.VF.height = height
import trimesh, time

def run():
    start = time.perf_counter()
    r = read(ROOT / 'FROZEN_COHORT.json')['rows'][1]
    a = npz(V4 / 'payload/whole_private' / r['key'] / 'reference.npz')
    p = npz(V4 / 'payload/whole_inputs' / r['key'] / 'preparation.npz')
    tri = a['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
    v = tri.reshape(-1, 3)
    f = np.arange(len(v)).reshape(-1, 3)
    roles = np.zeros(len(f), int)
    rec = dict(r, participant='REFERENCE_SELF_CONTROL')
    clean = inherited.VF.score_mesh(rec, v, f, roles)
    bad = inherited.VF.score_mesh(rec, v + np.array([5.0, 0, 0]), f, roles)
    assert clean['reconstruction_p95_mm'] < 1e-09 and bad['reconstruction_p95_mm'] > 0.35
    walls = []
    for width in [2.0, 3.5]:
        outer = trimesh.creation.box(extents=[4.0, 4.0, 4.0])
        inner = trimesh.creation.box(extents=[width] * 3)
        vv = np.r_[outer.vertices, inner.vertices]
        ff = np.r_[outer.faces, inner.faces + len(outer.vertices)]
        rr = np.r_[np.zeros(len(outer.faces), int), np.ones(len(inner.faces), int)]
        w = inherited.wall(vv, ff, rr)
        walls.append(dict(inner_width_mm=width, wall=w, passes_sampled=w['sampled_minimum_mm'] >= 0.5))
    assert walls[0]['passes_sampled'] and (not walls[1]['passes_sampled'])
    first = next((x for x in read(DATA / 'A/RECORDS.json') if x['status'] == 'EXPORTED'))
    m = npz(DATA / 'A' / first['participant'] / first['key'] / 'mesh.npz')
    pr = read(ROOT / 'PREREG_A.json')
    clean_row = inherited.score_row(first, m, pr)
    broken = dict(m)
    broken['faces'] = m['faces'][1:]
    broken['face_roles'] = m['face_roles'][1:]
    open_row = inherited.score_row(first, broken, pr)
    shifted = dict(m)
    shifted['vertices'] = m['vertices'] + np.array([30.0, 0.0, 1.0])
    moved_row = inherited.score_row(first, shifted, pr)
    assert clean_row['gates']['closed_single_component'] and (not open_row['gates']['closed_single_component'])
    assert clean_row['gates']['margin'] and (not moved_row['gates']['margin'])
    assert not moved_row['gates']['mesial_proximity'] and (not moved_row['gates']['distal_proximity'])
    mesh_faults = dict(source_key=first['key'], clean_gates=clean_row['gates'], missing_triangle_gates=open_row['gates'], translated_30x_1z_gates=moved_row['gates'], translation_mm=[30.0, 0.0, 1.0], scope='actual full inherited scorer, all original thresholds')
    result = dict(anatomy=dict(native_self_p95_mm=clean['reconstruction_p95_mm'], translated_5mm_p95_mm=bad['reconstruction_p95_mm'], fault_rejected=True, source=str(V4 / 'payload/whole_private' / r['key'] / 'reference.npz'), scope='anatomy-only selftest; open native scan not a closed-shell pass'), wall=walls, mesh_faults=mesh_faults, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/SCORER_CONTROLS.json', result)
    return result
if __name__ == '__main__':
    run()
