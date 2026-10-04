import os, sys, json, time, resource
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy import ndimage
import trimesh
sys.path.insert(0, '/oldcode')
sys.path.insert(0, '/oldcode/vendor')
sys.path.insert(0, '/deps')
from generate_r2 import section_mesh, largest_field, export, plain, dump
from generate import mesh2sdf, marching_cubes, height
import shapely
from shapely.geometry import Polygon

def close_outer(p):
    t = p['prior_vertices'][p['prior_faces'][p['prior_roles'] == 0]]
    m = trimesh.Trimesh(t.reshape(-1, 3), np.arange(len(t) * 3).reshape(-1, 3), process=True)
    from generate_r2 import loops
    section = Polygon()
    skirt = []
    for lp in loops(m):
        if abs(m.vertices[lp, 2] - float(p['margin_z'])).max() > 1e-06:
            raise ValueError('nonplanar inherited exterior boundary')
        section = section.symmetric_difference(Polygon(m.vertices[lp, :2]))
        for (i, j) in zip(lp, np.roll(lp, -1)):
            (a, b) = m.vertices[[i, j]]
            aa = a.copy()
            bb = b.copy()
            aa[2] -= 1.5
            bb[2] -= 1.5
            skirt.extend([[b, a, aa], [b, aa, bb]])
    base = []
    for t in shapely.constrained_delaunay_triangles(section).geoms:
        q = np.c_[np.asarray(t.exterior.coords)[:3], np.full(3, float(p['margin_z']) - 1.5)]
        if np.cross(q[1] - q[0], q[2] - q[0])[2] > 0:
            q = q[::-1]
        base.append(q)
    tri = np.r_[m.triangles, np.array(skirt), np.array(base)]
    m = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
    trimesh.repair.fix_normals(m)
    if not m.is_watertight:
        raise ValueError('outer cap failed')
    return m

def deform(p, m, method):
    xy = p['contact_xy']
    ceiling = p['contact_ceiling_mm']
    g = p['contact_gap_mm']
    roof = height(m.triangles, xy)
    valid = np.isfinite(g) & np.isfinite(roof) & np.isfinite(ceiling)
    near = valid & ((g <= 0.15) | (ceiling - roof <= 0.15))
    if not near.any():
        return (m, dict(active_nodes=0, max_displacement_mm=0.0))
    delta = np.clip(ceiling - g - roof, -4, 4)
    if method == 'scalar_control':
        delta[valid] = np.mean(delta[near])
    d = cKDTree(xy[near]).query(xy, workers=1)[0]
    blend = np.where(d < 1.0, 0.5 * (1 + np.cos(np.pi * np.minimum(d, 1))), 0.0)
    delta = np.where(valid, delta * blend, 0.0)
    interp = LinearNDInterpolator(xy, delta, fill_value=0.0)
    vv = m.vertices.copy()
    up = np.asarray(interp(vv[:, :2]))
    top = height(m.triangles, vv[:, :2])
    denom = top - float(p['margin_z'])
    w = np.clip((vv[:, 2] - float(p['margin_z'])) / np.maximum(denom, 0.1), 0, 1)
    vv[:, 2] += w * up
    derivative = 1 + up / np.maximum(denom, 0.1)
    if derivative.min() <= 0:
        raise ValueError('nonmonotone column deformation ' + str(derivative.min()))
    m.vertices = vv
    return (m, dict(active_nodes=int(near.sum()), max_displacement_mm=float(abs(up * w).max()), column_derivative_min=float(derivative.min()), conditioning='Added measured bounded gap sheet', global_self_intersection_proof='MISSING'))

def shell_from_outer(p, mesh):
    low = mesh.bounds[0] - 0.6
    high = mesh.bounds[1] + 0.6
    low[2] -= 0.6
    center = (low + high) / 2
    side = float((high - low).max())
    n = int(np.ceil(side / 0.15))
    h = side / n
    if n > 256:
        raise ValueError('grid budget exceeded')
    origin = center - side / 2
    axes = [origin[i] + h * np.arange(n) for i in range(3)]
    outer = mesh2sdf.compute((mesh.vertices - center) / (side / 2), mesh.faces, size=n, fix=False) * (side / 2)
    (inside, component) = largest_field(outer < 0)
    outer[~inside] = np.maximum(outer[~inside], h)
    erosion = 0.85 - ndimage.distance_transform_edt(inside) * h
    iz = max(1, int(np.floor((float(p['margin_z']) - origin[2]) / h)))
    monotone = erosion.copy()
    monotone[:, :, iz:] = np.maximum.accumulate(erosion[:, :, iz:], axis=2)
    monotone[:, :, :iz] = monotone[:, :, iz, None]
    (mask, icomp) = largest_field(monotone < 0)
    monotone[~mask] = np.maximum(monotone[~mask], h)
    monotone[:, :, 0] = h
    cavity_volume = float(((monotone < 0) & (axes[2][None, None, :] >= float(p['margin_z']))).sum() * h ** 3)
    if cavity_volume < 1:
        raise ValueError('cavity below1mm3')
    (eo, poly, oi) = section_mesh(outer, origin, h, float(p['margin_z']))
    (ei, ipoly, ii) = section_mesh(monotone, origin, h, float(p['margin_z']))
    if not poly.contains(ipoly):
        raise ValueError('inner section not contained')
    base = []
    for t in shapely.constrained_delaunay_triangles(poly.difference(ipoly)).geoms:
        q = np.c_[np.asarray(t.exterior.coords)[:3], np.full(3, float(p['margin_z']))]
        if np.cross(q[1] - q[0], q[2] - q[0])[2] > 0:
            q = q[::-1]
        base.append(q)
    tri = np.r_[eo.triangles, ei.triangles[:, ::-1], np.array(base)]
    roles = np.r_[np.zeros(len(eo.faces), int), np.ones(len(ei.faces), int), np.full(len(base), 2, int)]
    final = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)
    final.merge_vertices(digits_vertex=8)
    trimesh.repair.fix_normals(final)
    if not final.is_watertight or not final.is_winding_consistent or final.volume <= 0 or (len(final.split(only_watertight=False)) != 1):
        raise ValueError('joint shell topology failed')
    return (final.vertices, final.faces, roles, dict(grid_mm=h, grid_n=n, cavity_volume_mm3=cavity_volume, information_track='MEASURED_CONTACT_JOINT_VIRTUAL_PREPARATION', physical_fit='UNKNOWN', wall_status='requires final geometry validation'))

def main():
    out = Path('/output')
    records = json.loads(Path('/inputs/RECORDS.json').read_text())
    rows = []
    start = time.perf_counter()
    for r in records:
        with np.load(Path('/inputs') / (r['key'] + '.npz')) as z:
            p = {k: z[k] for k in z.files if k != 'contact_reference_gap_unclipped'}
        for method in ['spatial_sheet', 'scalar_control']:
            row = dict(key=r['key'], case_key=r['case_key'], family=r['family'], dataset=r['dataset'], source_fdi=r['source_fdi'], participant=method, status='FAILED')
            t = time.perf_counter()
            try:
                outer = close_outer(p)
                (outer, di) = deform(p, outer, method)
                (v, f, roles, gi) = shell_from_outer(p, outer)
                export(out / method / r['key'], v, f, roles, p, dict(deformation=di, **gi))
                row.update(status='EXPORTED', build=gi, deformation=di)
            except Exception as e:
                row['reason'] = repr(e)
            row['seconds'] = time.perf_counter() - t
            rows.append(row)
            dump(out / 'RECORDS.json', rows)
            print(row['key'], method, row['status'], round(row['seconds'], 2), row.get('reason', ''), flush=True)
    dump(out / 'COST.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    main()
