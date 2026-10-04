"""One local mesh -> signed-distance grid -> CSG graft -> STL export example."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import argparse, json, resource, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage
from skimage.measure import marching_cubes
from freeze import ROOT, DATA, sha, write
from planner import graft_levelset, stl

def run(cid='001', round_id='R2'):
    tick = time.perf_counter()
    plan = json.loads((DATA / round_id / cid / 'plan.json').read_text())
    if plan['status'] != 'PREDICTED':
        raise ValueError('illustrative case has no plan')
    with np.load(plan['arrays']['path']) as z:
        nodes = z['dp_nodes'].astype(float)
    radius = plan['grafts']['dp']['radius_mm']
    cache = Path(_release_expand('@DENTAL_WORK_ROOT@/G_mandible_postop/R2/components')) / f'{cid}_Pre.npz'
    with np.load(cache) as z:
        v = z['jaw_vertices'].astype(float)
        f = z['jaw_faces'].astype('i4')
    minimum = np.minimum(v.min(0), nodes.min(0) - radius) - 6
    maximum = np.maximum(v.max(0), nodes.max(0) + radius) + 6
    pitch = 1.5
    shape = np.ceil((maximum - minimum) / pitch).astype(int) + 1
    while np.prod(shape) > 800000:
        pitch *= 1.1
        shape = np.ceil((maximum - minimum) / pitch).astype(int) + 1
    enginefile = Path(os.environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'references/field_engine/src/field_engine/faltkarna_v1_mesh_to_sdf.py'
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(enginefile.parent))
    import faltkarna_v1_mesh_to_sdf as engine
    (inside, diag) = engine.solid_via_stralvindning(v, f, pitch, minimum, tuple(shape), cap=1000000, chunk_vox=100000)
    phi_pre = ndimage.distance_transform_edt(~inside, sampling=pitch) - ndimage.distance_transform_edt(inside, sampling=pitch)
    phi_pre -= np.sign(phi_pre) * pitch / 2
    indices = np.indices(tuple(shape), dtype='f4').reshape(3, -1).T
    points = minimum + pitch * indices
    phi_graft = graft_levelset(points, nodes, radius).reshape(tuple(shape))
    combined = np.minimum(phi_pre, phi_graft)
    (vertices, faces, _, _) = marching_cubes(combined.astype('f4'), 0, spacing=(pitch, pitch, pitch), allow_degenerate=False)
    vertices += minimum
    out = DATA / ('field_demo' if round_id == 'R1' else 'field_demo_' + round_id)
    out.mkdir(exist_ok=True)
    np.savez_compressed(out / 'fields.npz', pre=phi_pre.astype('f4'), graft=phi_graft.astype('f4'), union=combined.astype('f4'), origin=minimum, pitch=pitch)
    stl(out / 'Pre_plus_graft_surrogate.stl', (vertices, faces))
    record = dict(case=cid, round=round_id, mesh_to_sdf='native field engine ray winding -> EDT grid; approximate signed distance', graft_field='analytic cylinder/halfspace sign field; not Euclidean distance at cut edges', union='min(phi_Pre,phi_graft)', radius_mm=radius, grid_pitch_mm=pitch, shape=shape.tolist(), native_winding=diag, pre_inside_voxels=int(inside.sum()), pre_inside_voxels_lost=int(np.sum(inside & (combined >= 0))), engine_path=str(enginefile), engine_sha256=sha(enginefile), plan_sha256=sha(DATA / round_id / cid / 'plan.json'), exports=[dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()], wall_seconds=time.perf_counter() - tick, peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, manufacture_status='research polymer prototype only; not validated donor geometry or patient-specific cutting guide')
    if (ROOT / 'FIELD_DEMO.json').exists() and (not (ROOT / 'FIELD_DEMO_R1.json').exists()):
        (ROOT / 'FIELD_DEMO_R1.json').write_bytes((ROOT / 'FIELD_DEMO.json').read_bytes())
    write(ROOT / f'FIELD_DEMO_{round_id}.json', record)
    write(ROOT / 'FIELD_DEMO.json', record)
    return record
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='001')
    ap.add_argument('--round', default='R2')
    a = ap.parse_args()
    print(json.dumps(run(a.case, a.round), indent=2))
