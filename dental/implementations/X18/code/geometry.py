"""Fixed RAS geometry, cervical site and independent donor surfaces (mm)."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = str(min(4, int(os.environ.get(key, '1'))))
import hashlib, json, zipfile, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from continuous import broad_phase, extremum
H = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X18_crown_antagonist'))
ZIP = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'))
LABEL = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bits2Bites'))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False, default=lambda a: a.tolist() if isinstance(a, np.ndarray) else a.item() if isinstance(a, np.generic) else str(a)) + '\n')

def state(stage, gate, next_op, **kw):
    import datetime
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X18-crown-antagonist', status=stage, latest_gate=gate, next_operation=next_op, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), **kw))

def pair(case):
    out = {}
    man = {}
    with zipfile.ZipFile(ZIP) as z:
        for jaw in ['upper', 'lower']:
            name = next((n for n in z.namelist() if n.endswith('/' + str(case) + '/' + jaw + '.stl')))
            blob = z.read(name)
            count = int.from_bytes(blob[80:84], 'little')
            dt = np.dtype([('normal', '<f4', (3,)), ('v', '<f4', (3, 3)), ('attr', '<u2')])
            assert len(blob) == 84 + 50 * count
            tris = np.frombuffer(blob, dt, offset=84, count=count)['v'].astype(float)
            (v, inv) = np.unique(tris.reshape(-1, 3), axis=0, return_inverse=True)
            f = inv.reshape(-1, 3)
            lp = LABEL / str(case) / (jaw + '_labels.npz')
            meta = json.loads((LABEL / str(case) / 'labels+landmarks.json').read_text())['arches'][jaw]
            assert hashlib.sha256(blob).hexdigest() == meta['input_sha256']
            assert sha(lp) == meta['labels']['sha256']
            labs = np.load(lp)['labels']
            assert len(labs) == len(v)
            fl = labs[f]
            owner = np.where(fl[:, 0] == fl[:, 1], fl[:, 0], np.where(fl[:, 0] == fl[:, 2], fl[:, 0], np.where(fl[:, 1] == fl[:, 2], fl[:, 1], 0)))
            out[jaw] = dict(v=v, f=f, tri=tris, labels=labs, owner=owner, landmarks=meta['landmarks'])
            man[jaw] = dict(member=name, source_sha256=meta['input_sha256'], label_file=str(lp), label_sha256=sha(lp), source_bytes=len(blob))
    return (out, man)

def query_height(tri, xy, upper):
    """Exact projected affine height at vertices, radius-screened all source faces."""
    det = np.cross(tri[:, 1, :2] - tri[:, 0, :2], tri[:, 2, :2] - tri[:, 0, :2])
    ids = np.flatnonzero(np.abs(det) > 1e-12)
    t = tri[ids]
    center = t[:, :, :2].mean(1)
    radius = np.linalg.norm(t[:, :, :2] - center[:, None, :], axis=2).max(1)
    tree = cKDTree(center)
    hits = tree.query_ball_point(xy, float(radius.max()), workers=4)
    zz = np.full(len(xy), np.inf if upper else -np.inf)
    faces = np.full(len(xy), -1, int)
    for (k, js) in enumerate(hits):
        js = np.asarray(js, int)
        if not len(js):
            continue
        tt = t[js]
        (a, b, c) = (tt[:, 0, :2], tt[:, 1, :2], tt[:, 2, :2])
        den = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
        p = xy[k]
        w1 = ((b[:, 1] - c[:, 1]) * (p[0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (p[1] - c[:, 1])) / den
        w2 = ((c[:, 1] - a[:, 1]) * (p[0] - c[:, 0]) + (a[:, 0] - c[:, 0]) * (p[1] - c[:, 1])) / den
        ok = (w1 >= -1e-10) & (w2 >= -1e-10) & (1 - w1 - w2 >= -1e-10)
        if not ok.any():
            continue
        vals = w1[ok] * tt[ok, 0, 2] + w2[ok] * tt[ok, 1, 2] + (1 - w1 - w2)[ok] * tt[ok, 2, 2]
        ii = int(np.argmin(vals) if upper else np.argmax(vals))
        zz[k] = vals[ii]
        faces[k] = ids[js[ok][ii]]
    return (zz, faces)

def site(arch, fdi, n):
    """Only the cervical bottom30% is allowed to define target site."""
    pp = arch['v'][arch['labels'] == fdi]
    if len(pp) < 300:
        raise ValueError('Insufficient predicted tooth vertices')
    q = float(np.partition(pp[:, 2], int(0.3 * (len(pp) - 1)))[int(0.3 * (len(pp) - 1))])
    band = pp[pp[:, 2] <= q]
    (lo, hi) = np.quantile(band[:, :2], [0.02, 0.98], axis=0)
    center = (lo + hi) / 2
    half = (hi - lo) / 2
    if np.min(half) < 1.2 or np.max(half) > 8:
        raise ValueError('Invalid cervical site scale')
    ax = np.linspace(-1, 1, n)
    (x, y) = np.meshgrid(ax, ax, indexing='ij')
    uv = np.c_[x.ravel(), y.ravel()]
    keep = np.sum(uv ** 2, axis=1) <= 1 + 1e-12
    uv = uv[keep]
    xy = center + uv * half
    index = np.full((n, n), -1, int)
    index.ravel()[keep] = np.arange(keep.sum())
    faces = []
    for i in range(n - 1):
        for j in range(n - 1):
            for ij in [[(i, j), (i + 1, j), (i, j + 1)], [(i + 1, j), (i + 1, j + 1), (i, j + 1)]]:
                ids = [index[a, b] for (a, b) in ij]
                if min(ids) >= 0:
                    faces.append(ids)
    edges = np.concatenate([np.asarray(faces)[:, [0, 1]], np.asarray(faces)[:, [1, 2]], np.asarray(faces)[:, [2, 0]]])
    edges.sort(1)
    (eu, ec) = np.unique(edges, axis=0, return_counts=True)
    border = np.unique(eu[ec == 1])
    neighbor = [fdi - 1, fdi + 1]
    heights = []
    for f in neighbor:
        p = arch['v'][arch['labels'] == f]
        if len(p) >= 300:
            heights.append(float(np.quantile(p[:, 2], 0.95)))
    if not heights:
        raise ValueError('Missing usable adjacent tooth')
    return dict(center=center, half=half, z_cervical=q, z_top=float(np.mean(heights)), xy=xy, uv=uv, faces=np.asarray(faces, int), border=border, index=index, area_weight=np.full(len(xy), float(np.prod(2 * half / (n - 1)))), neighbors=neighbor)

def donor_template(pr):
    n = pr['grid_n']
    by = {p: [] for p in [4, 5, 6, 7]}
    manifest = {}
    start = time.perf_counter()
    for case in pr['donor_cases']:
        (pairdata, man) = pair(case)
        manifest[str(case)] = man
        arch = pairdata['lower']
        for fdi in pr['tooth_fdi']:
            try:
                s = site(arch, fdi, n)
                tri = arch['tri'][arch['owner'] == fdi]
                (zz, _) = query_height(tri, s['xy'], False)
                pp = arch['v'][arch['labels'] == fdi]
                top = float(np.quantile(pp[:, 2], 0.95))
                zz = zz - top
                zz[~np.isfinite(zz)] = np.nan
                by[fdi % 10].append(zz)
            except ValueError:
                continue
    template = {str(p): np.nanmedian(np.asarray(rows), axis=0) for (p, rows) in by.items()}
    for (p, t) in template.items():
        ok = np.isfinite(t)
        if not ok.all():
            t[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), t[ok])
    np.savez_compressed(D / 'donor_template.npz', **template)
    dump(H / 'raw/DONOR_FIT.json', dict(manifest=manifest, donors_per_type={str(p): len(a) for (p, a) in by.items()}, wall_seconds=time.perf_counter() - start, template_sha256=sha(D / 'donor_template.npz')))
    return template

def signed_gap(U, xy, z, faces, bounded=True):
    C = np.c_[xy, z][faces]
    pairs = broad_phase(U, C)
    if not len(pairs):
        return dict(minimum_gap_mm=None, maximum_penetration_mm=0, query_s=0)
    return extremum(U, C, pairs, bounded)

def repair(U, s, z, clearance):
    info = signed_gap(U, s['xy'], z, s['faces'])
    shift = max(0, clearance - (info['minimum_gap_mm'] or clearance))
    return (z - shift, dict(before=info, global_downward_repair_mm=shift, after=signed_gap(U, s['xy'], z - shift, s['faces'])))

def write_stl(p, vertices, faces):
    tri = vertices[faces]
    norm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm /= np.maximum(np.linalg.norm(norm, axis=1)[:, None], 1e-15)
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    a = np.zeros(len(faces), dt)
    a['n'] = norm
    a['v'] = tri
    Path(p).write_bytes(b'X18 conditional crown roof; mm'.ljust(80, b' ') + len(faces).to_bytes(4, 'little') + a.tobytes())

def shell(xy, z, faces, thickness):
    nv = len(xy)
    v = np.r_[np.c_[xy, z], np.c_[xy, z - thickness]]
    f = np.r_[faces, faces[:, ::-1] + nv]
    oriented = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    sort = np.sort(oriented, axis=1)
    (un, ct) = np.unique(sort, axis=0, return_counts=True)
    boundary = {tuple(e) for e in un[ct == 1]}
    sides = []
    for (a, b) in oriented:
        if tuple(sorted([a, b])) in boundary:
            sides.extend([[a, b, b + nv], [a, b + nv, a + nv]])
    return (v, np.r_[f, np.asarray(sides)])
