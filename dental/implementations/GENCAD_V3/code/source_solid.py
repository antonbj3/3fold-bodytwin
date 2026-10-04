"""R5: measured triangle patch, explicit basal closure, volume, prep, crown."""
import sys, time, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'vendor'))
import numpy as np, trimesh
from scipy import ndimage
from skimage.measure import marching_cubes
from util import *
from legacy.geometry import export_stl, boundary_edges
from preparation import bad_edges

def cap_patch(v, f, apex_z):
    e = boundary_edges(f)
    out = {}
    incoming = {}
    all_edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    if np.any(np.unique(all_edges, axis=0, return_counts=True)[1] > 2):
        raise ValueError('nonmanifold source patch')
    for (a, b) in e:
        a = int(a)
        b = int(b)
        if a in out or b in incoming:
            raise ValueError('branched source boundary')
        out[a] = b
        incoming[b] = a
    if set(out) != set(incoming):
        raise ValueError('source boundary not cycles')
    verts = v.tolist()
    faces = f.tolist()
    loops = []
    remaining = set(out)
    while remaining:
        start = min(remaining)
        loop = []
        cur = start
        while cur in remaining:
            loop.append(cur)
            remaining.remove(cur)
            cur = out[cur]
        if cur != start or len(loop) < 3:
            raise ValueError('invalid boundary loop')
        point = v[loop].mean(0)
        point[2] = apex_z
        index = len(verts)
        verts.append(point.tolist())
        for (a, b) in zip(loop, loop[1:] + loop[:1]):
            faces.append([b, a, index])
        loops.append(loop)
    vv = np.asarray(verts)
    ff = np.asarray(faces)
    if bad_edges(ff):
        raise ValueError('capped source is not edge-closed')
    return (vv, ff, loops)

def build():
    import crown_fit_geometry as G
    spec = importlib.util.spec_from_file_location('gencad_x1b_r5', ROOT / 'vendor/crown_design_geometry.py')
    CDG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CDG)
    CDG.G = G
    pr = read(ROOT / 'PREREG_R5.json')
    q = pr['parameters']
    prior = read(ROOT / 'FROZEN_PREPARATIONS_R3.json')['payload']
    meta = {r['case_key']: r for r in prior['records']}
    records = []
    for key in pr['cases']:
        start = time.perf_counter()
        record = dict(case_key=key, split=meta[key]['split'], status='FAILED_SOURCE_SOLID')
        try:
            source_path = PAYLOAD / 'preparations_r3' / key / 'fields.npz'
            if sha(source_path) != prior['files']['preparations_r3/' + key + '/fields.npz']['sha256']:
                raise ValueError('source drift')
            a = dict(np.load(source_path, allow_pickle=False))
            local = (a['source_vertices'] - a['center']) @ a['R']
            base = float(local[:, 2].min())
            (v, f, loops) = cap_patch(local, a['source_faces'], base - q['basal_apex_below_source_mm'])
            mesh = trimesh.Trimesh(v, f, process=False)
            h = q['grid_spacing_mm']
            grid = G.Grid(v.min(0) - q['padding_mm'], v.max(0) + q['padding_mm'], h)
            if np.prod(grid.shape) > q['max_grid_nodes']:
                raise ValueError('bounded grid node limit')
            occupied = np.zeros(grid.shape, bool)
            (xs, ys, zs) = grid.axes()
            (xx, yy) = np.meshgrid(xs, ys, indexing='ij')
            xy = np.c_[xx.ravel(), yy.ravel()]
            for (k, z) in enumerate(zs):
                pts = np.c_[xy, np.full(len(xy), z)]
                vals = []
                for block in np.array_split(pts, max(1, int(np.ceil(len(pts) / 2048)))):
                    vals.extend(mesh.contains(block))
                occupied[:, :, k] = np.asarray(vals).reshape(xx.shape)
            if not occupied.any():
                raise ValueError('empty capped-mesh volume')
            phi = (ndimage.distance_transform_edt(~occupied, sampling=h) - ndimage.distance_transform_edt(occupied, sampling=h)).astype(np.float32)
            margin = base + q['virtual_margin_above_min_source_mm']
            case = dict(grid=grid, phi_T=phi, z_m=margin)
            dest = PAYLOAD / 'preparations_r5' / key
            dest.mkdir(parents=True, exist_ok=True)
            scratch = dest / 'construction'
            scratch.mkdir(exist_ok=True)
            CDG.SCRATCH = str(scratch)
            (prep, mp) = CDG.preparation(case, q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
            (pv, pf) = G.surface_mesh(prep, grid)
            if np.einsum('ij,ij->i', pv[pf[:, 0]], np.cross(pv[pf[:, 1]], pv[pf[:, 2]])).sum() < 0:
                pf = pf[:, ::-1]
            export_stl(dest / 'preparation.stl', pv @ a['R'].T + a['center'], pf)
            export_stl(dest / 'native_exterior.stl', a['source_vertices'], a['source_faces'])
            np.savez_compressed(dest / 'fields.npz', source_vertices=a['source_vertices'], source_faces=a['source_faces'], phi=phi, preparation=prep, origin=grid.origin, step=h, R=a['R'], center=a['center'], preparation_vertices=pv, preparation_faces=pf, closed_source_vertices=v, closed_source_faces=f)
            for p in scratch.glob('*.npz'):
                p.unlink()
            scratch.rmdir()
            cavity = G.design_cavity(prep, grid, np.full(prep.shape, q['spacer_mm']), levels=(0.0,))
            z = grid.axes()[2][None, None, :]
            field = np.maximum.reduce([phi, -cavity, np.broadcast_to(margin - z, phi.shape)]).astype(np.float32)
            (cv, cf, _, _) = marching_cubes(field, q['iso_level_mm'], spacing=(h,) * 3, allow_degenerate=False)
            cv += grid.origin
            if np.einsum('ij,ij->i', cv[cf[:, 0]], np.cross(cv[cf[:, 1]], cv[cf[:, 2]])).sum() < 0:
                cf = cf[:, ::-1]
            cd = PAYLOAD / 'crowns_r5' / key
            cd.mkdir(parents=True, exist_ok=True)
            export_stl(cd / 'crown.stl', cv @ a['R'].T + a['center'], cf)
            np.savez_compressed(cd / 'prediction.npz', field=field, cavity=cavity, vertices=cv, faces=cf, margin_z=margin)
            record.update(status='EXPORTED', boundary_loops=len(loops), source_faces=len(a['source_faces']), closure_faces=len(f) - len(a['source_faces']), grid_shape=grid.shape, virtual_margin_z_mm=margin)
        except (ValueError, RuntimeError, IndexError) as e:
            record.update(reason=str(e))
        record['seconds'] = time.perf_counter() - start
        records.append(record)
        print('source-solid', key, record['status'], record.get('reason', ''), round(record['seconds'], 2), flush=True)
        budget()
    files = {str(p.relative_to(PAYLOAD)): {'sha256': sha(p), 'bytes': p.stat().st_size} for d in ['preparations_r5', 'crowns_r5'] for p in sorted((PAYLOAD / d).rglob('*')) if p.is_file()}
    freeze(ROOT / 'FROZEN_CROWNS_R5.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R5.json'), files=files, records=records, measurement_status='Source measurement held until frozen construction'))

class Remap:

    def __init__(self, b):
        self.base = b
        self.payload = b.payload / 'round5'

    def name(self, n):
        return n.replace('PREREG_R4.json', 'PREREG_R5.json').replace('FROZEN_CROWNS.json', 'FROZEN_CROWNS_R5.json').replace('payload/preparations_r3/', 'payload/preparations_r5/').replace('payload/crowns/', 'payload/crowns_r5/')

    def json(self, n):
        return self.base.json(self.name(n))

    def bytes(self, n):
        return self.base.bytes(self.name(n))

    def npz(self, n):
        return self.base.npz(self.name(n))

def evaluate(bundle):
    from crown_track import evaluate as check
    return check(Remap(bundle))
if __name__ == '__main__':
    build()
