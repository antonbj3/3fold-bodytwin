import thread_guard
from cadlib import *
from readonly import parents, scope
from pipeline import stitch
from margin import cycles
import trimesh, three_mf, resource
P = parents()
with scope(BASE / 'PROOF_LANE_FULL_CROWN_R4/code', ['common', 'geometry', 'generate', 'generate_c', 'generate_d']):
    AD = module('generate_d', BASE / 'PROOF_LANE_FULL_CROWN_R4/code/generate_d.py')

def predict():
    st = time.perf_counter()
    inputs = {r['key']: r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']}
    profiles = {r['id']: r for r in read(ROOT / 'MATERIAL_PROFILES.json')['profiles']}
    rows = []
    for row in read(ROOT / 'raw/R4B_PREDICTIONS.json'):
        r = dict(row)
        if r['status'] != 'GENERATED':
            rows.append(r)
            continue
        a = np.load(r['mesh_path'])
        ext = trimesh.Trimesh(a['ov'], a['of'], process=False)
        n = np.load(inputs[r['key']]['public_geometry'])['neighbors']
        mid = n.mean(1)[:, 0]
        (e, info) = AD.adapt(ext, dict(mesial=n[mid < 0], distal=n[mid >= 0]))
        base = e.vertices.copy()
        z0 = float(a['iv'][:, 2].min())
        ob = np.flatnonzero(abs(base[:, 2] - z0) < 1e-07)
        ib = np.arange(24)
        attempts = []
        accepted = None
        for alpha in np.arange(9) * 0.25:
            ov = base.copy()
            ov[ob, :2] += alpha * ov[ob, :2] / np.linalg.norm(ov[ob, :2], axis=1)[:, None]
            w = P['cap'].wall_check(ov, e.faces, a['iv'], a['inf'], profiles[r['material']]['wall_mm'])
            attempts.append(dict(extra_flare_mm=float(alpha), wall=w))
            if w['status'] == 'PASS':
                accepted = (ov, float(alpha))
                break
        if accepted is None:
            ov = base.copy()
            alpha = 0.0
            repair_status = 'NO_WALL_FEASIBLE_IN_FROZEN_GRID'
        else:
            (ov, alpha) = accepted
            repair_status = 'WALL_FEASIBLE_DISCRETE_RING'
        (m, roles) = stitch(ov, e.faces, a['iv'], a['inf'], ib, ob)
        dest = DATA / 'R5' / r['uid']
        dest.mkdir(exist_ok=True, parents=True)
        npz = dest / 'crown.npz'
        np.savez_compressed(npz, vertices=m.vertices, faces=m.faces, roles=roles, ov=ov, of=e.faces, iv=a['iv'], inf=a['inf'], roof=a['roof'], prior=a['prior'], undercut_dot=a['undercut_dot'], undercut_faces=a['undercut_faces'])
        stl = dest / 'crown.stl'
        m.export(stl)
        mf = dest / 'crown.3mf'
        three_mf.write(mf, m.vertices, m.faces, roles, dict(case_id=r['key'], material=r['material'], physical_release='UNKNOWN', parent_sha256=r['mesh_sha256']))
        r.update(parent_mesh_sha256=r['mesh_sha256'], mesh_path=npz, mesh_sha256=sha(npz), stl_path=stl, stl_sha256=sha(stl), three_mf_path=mf, three_mf_sha256=sha(mf), neighbour_adaptation=info, ring_repair=dict(status=repair_status, extra_flare_mm=alpha, attempts=attempts))
        rows.append(r)
        dump(ROOT / 'raw/R5_PARTIAL.json', rows)
        print(r['uid'], repair_status, alpha, flush=True)
    dump(ROOT / 'raw/R5_PREDICTIONS.json', rows)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R5.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R5.json'), rows_sha256=sha(ROOT / 'raw/R5_PREDICTIONS.json'), code_sha256=sha(Path(__file__)), parent_prediction_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R4B.json'), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, external_reference_scored=False))
    state('R5_PREDICTIONS_FROZEN', dict(requested=len(rows)), 'Score coupled context after real neighbour deformation and wall repair')
if __name__ == '__main__':
    predict()
