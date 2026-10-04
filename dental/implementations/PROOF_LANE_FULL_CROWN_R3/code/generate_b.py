"""Preserve the measured exterior boundary; use the SDF only for its cavity."""
import resource
resource.setrlimit(resource.RLIMIT_AS, (3584 * 1024 ** 2, 3584 * 1024 ** 2))
import generate_a as ga
from memory_height import height
from generate_r2 import loops, section_mesh, largest_field
from generate import export, dump, mesh2sdf
from pathlib import Path
import sys, json, time
import numpy as np
import trimesh, shapely
from shapely.geometry import Polygon
from scipy import ndimage

def construct(p):
    xy = p['contact_xy']
    g = p['contact_gap_mm']
    ceiling = p['contact_ceiling_mm']
    margin = float(p['margin_z'])
    prior = ga.close_outer(p)
    z = height(prior.triangles, xy)
    measured = np.isfinite(g) & np.isfinite(ceiling)
    z[measured] = ceiling[measured] - g[measured]
    faces = p['contact_faces']
    keep = np.all(np.isfinite(z[faces]) & (z[faces] > margin), axis=1)
    faces = faces[keep]
    used = np.unique(faces)
    index = np.full(len(xy), -1, int)
    index[used] = np.arange(len(used))
    top = trimesh.Trimesh(np.c_[xy[used], z[used]], index[faces], process=False)
    parts = top.split(only_watertight=False)
    if len(parts) != 1:
        raise ValueError('top domain not single component:' + str(len(parts)))
    poly = Polygon()
    skirt = []
    for loop in loops(top):
        poly = poly.symmetric_difference(Polygon(top.vertices[loop, :2]))
        for (i, j) in zip(loop, np.roll(loop, -1)):
            (a, b) = top.vertices[[i, j]]
            aa = a.copy()
            bb = b.copy()
            aa[2] = margin
            bb[2] = margin
            skirt.extend([[b, a, aa], [b, aa, bb]])
    if not poly.is_valid:
        raise ValueError('invalid projected domain')
    exterior = np.r_[top.triangles, np.asarray(skirt)]
    m = trimesh.Trimesh(exterior.reshape(-1, 3), np.arange(len(exterior) * 3).reshape(-1, 3), process=True)
    skirtbelow = []
    for loop in loops(m):
        for (i, j) in zip(loop, np.roll(loop, -1)):
            (a, b) = m.vertices[[i, j]]
            aa = a.copy()
            bb = b.copy()
            aa[2] -= 1.5
            bb[2] -= 1.5
            skirtbelow.extend([[b, a, aa], [b, aa, bb]])
    bottom = []
    for t in shapely.constrained_delaunay_triangles(poly).geoms:
        q = np.c_[np.asarray(t.exterior.coords)[:3], np.full(3, margin - 1.5)]
        if np.cross(q[1] - q[0], q[2] - q[0])[2] > 0:
            q = q[::-1]
        bottom.append(q)
    tri = np.r_[exterior, np.array(skirtbelow), np.array(bottom)]
    solid = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
    trimesh.repair.fix_normals(solid)
    if not solid.is_watertight:
        raise ValueError('explicit outer solid not closed')
    low = solid.bounds[0] - 0.6
    high = solid.bounds[1] + 0.6
    low[2] -= 0.6
    center = (low + high) / 2
    side = float((high - low).max())
    n = int(np.ceil(side / 0.15))
    h = side / n
    if n > 256:
        raise ValueError('grid budget')
    origin = center - side / 2
    field = mesh2sdf.compute((solid.vertices - center) / (side / 2), solid.faces, size=n, fix=False) * (side / 2)
    (inside, _) = largest_field(field < 0)
    cavity = (0.85 - ndimage.distance_transform_edt(inside) * h).astype(np.float32)
    del field, inside
    iz = max(1, int(np.floor((margin - origin[2]) / h)))
    cavity[:, :, iz:] = np.maximum.accumulate(cavity[:, :, iz:], axis=2)
    cavity[:, :, :iz] = cavity[:, :, iz, None]
    (mask, component) = largest_field(cavity < 0)
    cavity[~mask] = np.maximum(cavity[~mask], h)
    cavity[:, :, 0] = h
    cv = float(((cavity < 0) & ((origin[2] + h * np.arange(n))[None, None, :] >= margin)).sum() * h ** 3)
    if cv < 1:
        raise ValueError('cavity below1mm3')
    (inner, ipoly, info) = section_mesh(cavity, origin, h, margin)
    del cavity, mask
    if not poly.contains(ipoly):
        raise ValueError('cavity outside explicit boundary')
    annulus = []
    for t in shapely.constrained_delaunay_triangles(poly.difference(ipoly)).geoms:
        q = np.c_[np.asarray(t.exterior.coords)[:3], np.full(3, margin)]
        if np.cross(q[1] - q[0], q[2] - q[0])[2] > 0:
            q = q[::-1]
        annulus.append(q)
    tri = np.r_[exterior, inner.triangles[:, ::-1], np.array(annulus)]
    roles = np.r_[np.zeros(len(exterior), int), np.ones(len(inner.faces), int), np.full(len(annulus), 2, int)]
    final = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)
    final.merge_vertices(digits_vertex=8)
    trimesh.repair.fix_normals(final)
    if not final.is_watertight or not final.is_winding_consistent or final.volume <= 0 or (len(final.split(only_watertight=False)) != 1):
        raise ValueError('explicit joint shell topology failed')
    pred = height(final.vertices[final.faces[roles == 0]], xy)
    good = measured & np.isfinite(pred)
    mismatch = np.abs(pred[good] - (ceiling - g)[good])
    return (final.vertices, final.faces, roles, dict(grid_n=n, grid_mm=h, cavity_volume_mm3=cv, represented_measured_nodes=int(good.sum()), measurement_nodes=int(measured.sum()), node_preservation_max_mm=float(mismatch.max()) if len(mismatch) else None, domain_dropped_faces=int((~keep).sum()), information_track='MEASURED_CONTACT_JOINT_VIRTUAL_PREPARATION', physical_fit='UNKNOWN', axial_surface='vertical skirt closure, not measured anatomy', global_self_intersection_proof='MISSING'))

def main():
    key = sys.argv[1]
    out = Path('/output')
    r = next((r for r in json.loads(Path('/inputs/RECORDS.json').read_text()) if r['key'] == key))
    row = {k: r[k] for k in ['key', 'case_key', 'family', 'dataset', 'source_fdi']}
    row.update(participant='explicit_sheet', status='FAILED')
    start = time.perf_counter()
    try:
        with np.load(Path('/inputs') / (key + '.npz')) as z:
            p = {k: z[k] for k in z.files if k != 'contact_reference_gap_unclipped'}
        (v, f, roles, info) = construct(p)
        export(out / 'explicit_sheet' / key, v, f, roles, p, info)
        row.update(status='EXPORTED', build=info)
    except Exception as e:
        row['reason'] = repr(e)
        import traceback
        traceback.print_exc()
    row.update(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(out / (key + '.json'), row)
    print(json.dumps(row), flush=True)
if __name__ == '__main__':
    main()
