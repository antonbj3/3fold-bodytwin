"""Whole measured surface -> complete research crown/die; separate full-prescan information track."""
from dental_release.paths import expand as _release_expand
from fc_common import *
import trimesh, shapely, resource
from shapely.geometry import Polygon, Point
sys.path.insert(0, _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'))
from homolog_plugin import main_component, cap

def down(x):
    return np.nextafter(x, -np.inf)

def up(x):
    return np.nextafter(x, np.inf)

def norm_interval(v):
    lo = down(v)
    hi = up(v)
    amin = np.where((lo <= 0) & (hi >= 0), 0, np.minimum(abs(lo), abs(hi)))
    amax = np.maximum(abs(lo), abs(hi))
    ls = down(amin * amin)
    hs = up(amax * amax)
    sl = down(down(ls[..., 0] + ls[..., 1]) + ls[..., 2])
    sh = up(up(hs[..., 0] + hs[..., 1]) + hs[..., 2])
    return (down(np.sqrt(np.maximum(sl, 0))), up(np.sqrt(np.maximum(sh, 0))))

def wall_bound(tri, c, r=0.0):
    q = tri - c
    q[:, :, 2] = np.maximum(q[:, :, 2], 0)
    (lv, _) = norm_interval(q)
    diam = np.zeros(len(tri))
    for (i, j) in [(0, 1), (1, 2), (2, 0)]:
        diam = np.maximum(diam, norm_interval(tri[:, i] - tri[:, j])[1])
    lower = down(down(lv.min(1) - diam) - r)
    return (float(lower.min()), lower)

def loops_of(m):
    f = m.faces
    ed = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    (_, ix, ct) = np.unique(np.sort(ed, axis=1), axis=0, return_inverse=True, return_counts=True)
    pending = {tuple(e) for e in ed[ct[ix] == 1].tolist()}
    loops = []
    while pending:
        (a, b) = pending.pop()
        loop = [a, b]
        while loop[-1] != loop[0]:
            e = next((e for e in pending if e[0] == loop[-1]), None)
            if e is None:
                raise ValueError('boundary not a directed loop')
            pending.remove(e)
            loop.append(e[1])
        loops.append(loop[:-1])
    return loops

def cavity_triangles(c, r, margin, closed=False):
    N = 128
    M = 32
    theta = 2 * np.pi * np.arange(N) / N
    ringxy = np.c_[np.cos(theta), np.sin(theta)]
    v = []
    for (z, rad) in [(margin, r)] + [(c[2] + r * np.sin(j * np.pi / (2 * M)), r * np.cos(j * np.pi / (2 * M))) for j in range(M)]:
        v.extend(np.c_[c[:2] + rad * ringxy, np.full(N, z)].tolist())
    top = len(v)
    v.append([c[0], c[1], c[2] + r])
    f = []
    for layer in range(M):
        for k in range(N):
            a = layer * N + k
            b = layer * N + (k + 1) % N
            cc = (layer + 1) * N + k
            d = (layer + 1) * N + (k + 1) % N
            f.extend([[a, b, cc], [b, d, cc]])
    for k in range(N):
        f.append([M * N + k, M * N + (k + 1) % N, top])
    if closed:
        base = len(v)
        v.append([c[0], c[1], margin])
        for k in range(N):
            f.append([base, (k + 1) % N, k])
    m = trimesh.Trimesh(np.array(v), np.array(f), process=False)
    trimesh.repair.fix_normals(m)
    return m

def build(tri, margin):
    (m, drop, ncomp) = main_component(tri)
    (closed, ci) = cap(m)
    if not closed.is_watertight:
        raise ValueError('source closure not watertight')
    ext = closed.slice_plane([0, 0, margin], [0, 0, 1], cap=False)
    ext.merge_vertices()
    loops = loops_of(ext)
    if len(loops) != 1:
        raise ValueError(f'multiple cervical boundaries:{len(loops)}')
    ring = ext.vertices[loops[0]].copy()
    if np.max(abs(ring[:, 2] - margin)) > 1e-07:
        raise ValueError('unclosed source holes above virtual margin')
    ring[:, 2] = margin
    ext.vertices[loops[0], 2] = margin
    poly = Polygon(ring[:, :2])
    if not poly.is_valid:
        raise ValueError('invalid cervical polygon')
    bounds = poly.bounds
    xy = []
    for fx in np.linspace(0.2, 0.8, 9):
        for fy in np.linspace(0.2, 0.8, 9):
            q = np.array([bounds[0] + fx * (bounds[2] - bounds[0]), bounds[1] + fy * (bounds[3] - bounds[1])])
            if poly.contains(Point(q)):
                xy.append(q)
    xy += [np.array(poly.centroid.coords[0]), np.array(poly.representative_point().coords[0])]
    top = margin + 0.5 * (ext.vertices[:, 2].max() - margin)
    best = None
    for p in xy:
        c = np.r_[p, top]
        (bound, _) = wall_bound(ext.triangles, c)
        r = min(2.0, bound - 0.5 - 2e-05)
        if best is None or r > best[0]:
            best = (r, c, bound)
    (r, c, bound) = best
    if r < 0.43:
        raise ValueError(f'no certified cavity radius: r={r:.9g}mm')
    inner = cavity_triangles(c, r, margin)
    iv = inner.vertices
    iff = inner.faces[:, ::-1]
    inner_ring = iv[:128, :2]
    annulus = Polygon(ring[:, :2], [inner_ring])
    if not annulus.is_valid or not poly.contains(Polygon(inner_ring)):
        raise ValueError('cavity rim outside observed section')
    triangles = shapely.constrained_delaunay_triangles(annulus)
    base = []
    for g in triangles.geoms:
        pp = np.array(g.exterior.coords)[:3]
        tt = np.c_[pp, np.full(3, margin)]
        if np.cross(tt[1] - tt[0], tt[2] - tt[0])[2] > 0:
            tt = tt[::-1]
        base.append(tt)
    tr = np.concatenate([ext.triangles, iv[iff], np.array(base)])
    roles = np.r_[np.zeros(len(ext.faces), int), np.ones(len(iff), int), np.full(len(base), 2, int)]
    mesh = trimesh.Trimesh(tr.reshape(-1, 3), np.arange(len(tr) * 3).reshape(-1, 3), process=False)
    mesh.merge_vertices()
    trimesh.repair.fix_normals(mesh)
    assert len(mesh.faces) == len(roles)
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume <= 0:
        raise ValueError('stitched shell fails topology')
    (actual_bound, _) = wall_bound(mesh.triangles[roles == 0], c, r)
    if actual_bound < 0.5:
        raise ValueError('actual outer wall bound below0.5')
    die = cavity_triangles(c, r - 0.08, margin - 1.0, closed=True)
    angle = 355 / (113 * 32)
    film_lower = 0.08 * (1 - angle * angle / 2)
    film_upper = 0.08
    return (mesh, roles, die, dict(cavity_radius_mm=r, die_radius_mm=r - 0.08, axis_top=c.tolist(), certified_wall_lower_mm=actual_bound, normal_film_interval_mm=[film_lower, film_upper], minimum_spherical_tool_radius_mm=0.25, ideal_insertion='PROVED_DOWNWARD_CLOSED_CAVITY', ideal_spherical_tool_access='PROVED_CONNECTED_DOWNWARD_BALL_SWEEP', full_CAM='UNKNOWN_HOLDER_ORIENTATION_MACHINE', physical_seated_film='UNKNOWN', source_rejected_area_fraction=drop, source_components=ncomp, source_cap=ci, source_preservation='clipped measured triangles plus explicitly virtual closures', shell_components=len(mesh.split(only_watertight=False)), wall_enclosure='all exterior triangles; outward-rounded vertex ray-distance minus triangle diameter; no affine sensitivity'))

def run():
    st = time.perf_counter()
    out = DATA / 'R4_predictions'
    out.mkdir(exist_ok=True)
    records = read(V4 / 'payload/whole_inputs/RECORDS.json')
    rows = []
    for rec in records:
        row = dict(rec, participant='prescan_joint_die_shell', kind='shell', status='FAILED', information_track='FULL_PREOPERATIVE_SCAN_JOINT_PREPARATION')
        try:
            p = npz(V4 / 'payload/whole_inputs' / rec['key'] / 'preparation.npz')
            a = npz(V4 / 'payload/whole_private' / rec['key'] / 'reference.npz')
            (m, roles, die, info) = build(a['source_triangles'], float(p['margin_z']))
            dest = out / row['participant'] / rec['key']
            dest.mkdir(exist_ok=True, parents=True)
            np.savez_compressed(dest / 'mesh.npz', vertices=m.vertices, faces=m.faces, face_roles=roles)
            world = trimesh.Trimesh(m.vertices @ p['source_R'].T + p['source_base'], m.faces, process=False)
            world.export(dest / 'crown.stl')
            trimesh.Trimesh(die.vertices @ p['source_R'].T + p['source_base'], die.faces, process=False).export(dest / 'die.stl')
            dump(dest / 'CERTIFICATE.json', info)
            row.update(status='EXPORTED', certificate=info)
        except Exception as e:
            row['reason'] = str(e)
        rows.append(row)
        dump(out / 'RECORDS.json', rows)
        print(rec['key'], row['status'], row.get('reason', round(row.get('certificate', {}).get('cavity_radius_mm', 0), 4)), flush=True)
    dump(out / 'COST.json', dict(seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    freeze(ROOT / 'FROZEN_PREDICTIONS_R4.json', dict(files={str(p.relative_to(out)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'code/prescan_design.py'), information_track='FULL_PREOPERATIVE_SCAN_JOINT_PREPARATION', reference_queries='Full target is explicitly observed input; hidden-target generalization NOT CLAIMED', physical_measurement='NOT_RUN'))
if __name__ == '__main__':
    run()
