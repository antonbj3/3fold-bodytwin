"""Public preparation + neighbour surfaces only; withheld target never serialized."""
from common import *
import trimesh, mesh2sdf
from scipy.ndimage import map_coordinates
from skimage.measure import marching_cubes

def run():
    records = [r for r in read(V4 / 'payload/whole_inputs/RECORDS.json') if r['status'] == 'PREPARED'][:3]
    out = DATA / 'external_inputs'
    out.mkdir(exist_ok=True)
    locked = []
    for r in records:
        key = r['key']
        p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
        ctx = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
        (v, f, _, _) = marching_cubes(p['preparation'], 0, spacing=(float(p['step']),) * 3)
        v += p['origin']
        neighbors = np.concatenate([ctx['mesial_triangles'], ctx['distal_triangles']])
        tri = np.concatenate([v[f], neighbors])
        verts = tri.reshape(-1, 3)
        faces = np.arange(len(verts)).reshape(-1, 3)
        center = (verts.min(0) + verts.max(0)) / 2
        scale = 0.8 / np.max(np.linalg.norm(verts - center, axis=1))
        nv = (verts - center) * scale
        sdf = mesh2sdf.compute(nv, faces, size=64).astype(np.float32)
        coords = -1 + np.arange(64) * (2 / 64)
        xx = np.stack(np.meshgrid(coords, coords, coords, indexing='ij'), -1)
        world = xx / scale + center
        ix = ((world - p['origin']) / float(p['step'])).reshape(-1, 3).T
        cavity = map_coordinates(p['cavity'], ix, order=1, mode='constant', cval=100).reshape((64,) * 3)
        low = v.min(0) - [2, 2, 0]
        high = v.max(0) + [2, 2, 5]
        low[2] = float(p['margin_z'])
        np.savez_compressed(out / (key + '.npz'), tsdf=sdf, center=center, scale=scale, cavity=cavity, roi_low=low, roi_high=high, margin_z=p['margin_z'], source_R=p['source_R'], source_base=p['source_base'], preparation_vertices=v, preparation_faces=f, neighbor_triangles=neighbors)
        check = np.max(abs(nv / scale + center - verts))
        assert check < 1e-09
        locked.append(dict(**r, path=str(out / (key + '.npz')), sha256=sha(out / (key + '.npz')), roundtrip_error_mm=check, mm_per_voxel=2 / 64 / scale, context='public virtual preparation + measured neighbours, no original target surface', orientation='native site z-up; matching author ODD orientation is unverified'))
    freeze(ROOT / 'EXTERNAL_INPUTS.json', dict(rows=locked, selection='first three prepared records, all from one patient-case; three tooth classes', normalization='mesh2sdf grid [-1,1) with 2/64 spacing; max radius0.8'))
if __name__ == '__main__':
    run()
