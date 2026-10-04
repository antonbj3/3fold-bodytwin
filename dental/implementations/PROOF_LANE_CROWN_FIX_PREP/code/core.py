from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, time, datetime, resource, subprocess, importlib.util
from pathlib import Path
for x in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'IGL_NUM_THREADS']:
    os.environ[x] = '1'
os.environ['MPLBACKEND'] = 'Agg'
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
R = Path(os.environ.get('DENT_FIX_RUN_ROOT', Path(__file__).resolve().parents[1]))
B = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts'))
D = Path(os.environ.get('DENT_FIX_DATA_ROOT', _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_CROWN_FIX_PREP')))
PARENT_D = Path(_release_expand('@DENTAL_WORK_ROOT@'))
sys.path.insert(0, str(PARENT_D / 'PROOF_LANE_FULL_CROWN_R5/deps'))
import numpy as np, trimesh, igl
from scipy.optimize import linprog
from fractions import Fraction as F
from threadpoolctl import threadpool_limits
threadpool_limits(1)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def read(p):
    return json.loads(Path(p).read_text())

def clean(x):
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    q.replace(p)

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for s in iter(lambda : f.read(1024 ** 2), b''):
            h.update(s)
    return h.hexdigest()

def lock(name, x):
    p = R / name
    if p.exists():
        raise ValueError('Already frozen ' + name)
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def load(p):
    with np.load(p, allow_pickle=False) as a:
        return dict(a)

def module(name, p):
    sp = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(sp)
    sys.modules[name] = m
    sp.loader.exec_module(m)
    return m

def state(phase, latest, nxt):
    dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-fix-prep', updated_utc=now(), phase=phase, latest_gate=latest, next_operation=nxt, review_state='PENDING_INDEPENDENT_REVIEW', data_root=D))

def meshwrite(p, v, f):
    with open(p, 'w') as h:
        h.write(f'{len(v)} {len(f)}\n')
        np.savetxt(h, v, fmt='%.17g')
        np.savetxt(h, f, fmt='%d')

def meshread(p):
    with open(p) as h:
        (nv, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=nv)
        f = np.loadtxt(h, dtype=int)
    return (v, f.reshape(-1, 3))

def compact(v, f):
    (ids, inv) = np.unique(f, return_inverse=True)
    return (v[ids].copy(), inv.reshape(-1, 3))

def triangles_mesh(t):
    (v, inv) = np.unique(t.reshape(-1, 3), axis=0, return_inverse=True)
    m = trimesh.Trimesh(v, inv.reshape(-1, 3), process=False)
    trimesh.repair.fix_normals(m, multibody=True)
    return (m.vertices.copy(), m.faces.copy())

class Distance:

    def __init__(self, v, f):
        self.v = np.ascontiguousarray(v)
        self.f = np.ascontiguousarray(f)
        self.tree = igl.AABB()
        self.tree.init(self.v, self.f)

    def query(self, p):
        (d, j, q) = self.tree.squared_distance(self.v, self.f, np.ascontiguousarray(p))
        return (np.sqrt(np.maximum(0, d)), j, q)

    def inside(self, p):
        return np.abs(igl.fast_winding_number(self.v, self.f, np.ascontiguousarray(p))) > 0.5

def sample(t, n=4096):
    areas = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    idx = np.searchsorted(np.cumsum(areas), (np.arange(n) + 0.5) * areas.sum() / n)
    return t[idx].mean(1)

def shape_score(a, b):
    (av, af) = triangles_mesh(a)
    (bv, bf) = triangles_mesh(b)
    da = Distance(bv, bf).query(sample(a))[0]
    db = Distance(av, af).query(sample(b))[0]
    return dict(p95_mm=float(max(np.quantile(da, 0.95), np.quantile(db, 0.95))), sampled_max_mm=float(max(da.max(), db.max())), scope='area-stratified facet-centroid sample; no source uncertainty or full Hausdorff enclosure')

def close(v, f, key):
    p = D / (key + '_open.mesh')
    out = D / (key + '_closed.mesh')
    meshwrite(p, v, f)
    z = subprocess.run([str(PARENT_D / 'PROOF_LANE_CROWN_ROOT_CAUSE/close_native'), str(p), str(out)], capture_output=True, text=True, timeout=60)
    if z.returncode:
        raise ValueError('CLOSURE_EXIT ' + str(z.returncode) + ' ' + z.stderr[:100])
    info = json.loads(z.stdout)
    (vv, ff) = meshread(out)
    mm = trimesh.Trimesh(vv, ff, process=False)
    if mm.volume < 0:
        ff = ff[:, ::-1]
        meshwrite(out, vv, ff)
    info.update(path=str(out), sha256=sha(out), original_count=len(f), watertight=bool(mm.is_watertight))
    dump(R / 'raw' / f'{key}_CLOSURE.json', info)
    if info['intersection_count'] or not info['closed']:
        raise ValueError('INVALID_CLOSED_SUPPORT ' + json.dumps(info))
    return (vv, ff, info)

def boolean(a, b, out):
    z = subprocess.run([str(PARENT_D / 'PROOF_LANE_CROWN_ROOT_CAUSE/boolean'), str(a), str(b), str(out)], capture_output=True, text=True, timeout=90)
    if z.returncode:
        raise ValueError('BOOLEAN_EXIT ' + str(z.returncode) + ' ' + z.stdout[:120] + ' ' + z.stderr[:100])
    return json.loads(z.stdout)

def intersections(v, f, key):
    p = D / (key + '_check.mesh')
    meshwrite(p, v, f)
    z = subprocess.run([str(PARENT_D / 'PROOF_LANE_FULL_CROWN_R6/intersections'), str(p)], capture_output=True, text=True, timeout=60)
    if z.returncode:
        raise ValueError('INTERSECTION_EXIT ' + str(z.returncode))
    return json.loads(z.stdout)
K = 0.12
RHO = 0.5
G = 0.05
SEC = np.sqrt(1 + K * K)
SIN = K / SEC

def offset(z, H, rho=0.5):
    z = np.asarray(z)
    hc = H - rho
    tang = hc + rho * SIN
    return np.where(z <= tang, rho * SEC - K * (z - hc), np.sqrt(np.maximum(0, rho * rho - (z - hc) ** 2)))

def body(c, a, b, H, M, rho=0.5, bottom_override=None, n=64):
    tang = H - rho + rho * SIN
    bt = b if bottom_override is None else bottom_override
    zs = np.r_[np.linspace(bt, tang, max(2, int(np.ceil((tang - bt) / 0.1)) + 1)), H - rho + rho * np.sin(np.linspace(np.arctan(K), np.pi / 2, 17)[1:])]
    rad = a + offset(zs, H, rho)
    theta = np.arange(n) * 2 * np.pi / n
    circle = np.c_[np.cos(theta), np.sin(theta)]
    vv = np.concatenate([np.c_[c + rr * circle, np.full(n, z)] for (rr, z) in zip(rad, zs)]) @ M
    bot = len(vv)
    top = bot + 1
    v = np.r_[vv, np.array([[*c, bt], [*c, H]]) @ M]
    fs = []
    for j in range(len(zs) - 1):
        for i in range(n):
            k = (i + 1) % n
            fs.extend([[j * n + i, j * n + k, (j + 1) * n + k], [j * n + i, (j + 1) * n + k, (j + 1) * n + i]])
    for i in range(n):
        k = (i + 1) % n
        fs.extend([[bot, k, i], [(len(zs) - 1) * n + i, (len(zs) - 1) * n + k, top]])
    return (v, np.array(fs, int), zs, rad)

def frame_candidates(key):
    rows = read(R / 'raw/CERVICAL_FRAME_DIAGNOSTIC.json')
    out = [np.array(next((r['M'] for r in rows if r['key'] == key))), np.eye(3)]
    old = read(B / 'PROOF_LANE_CROWN_ROOT_CAUSE/RESULTS_D.json')['rows']
    q = next((r for r in old if r['key'] == key))
    if q['status'] == 'GENERATED':
        out.append(np.array(q['chosen']['world_to_local']))
    return out

def query_ball_rows(v, zs, n=64, w=0.0):
    rows = []
    wq = F(str(w))
    scale = 10 ** 6
    for j in range(len(zs) - 1):
        pts = v[j * n:(j + 2) * n]
        cen = (pts[:n].mean(0) + pts[n:].mean(0)) / 2
        cc = [F(float(x)) for x in cen]
        sq = max((sum(((F(float(x)) - c) ** 2 for (x, c) in zip(p, cc))) for p in pts))
        m = __import__('math').isqrt(sq.numerator * scale * scale // sq.denominator)
        while F(m * m, scale * scale) < sq:
            m += 1
        rr = F(m, scale) + wq
        r = float(rr)
        if F(r) < rr:
            r = np.nextafter(r, np.inf)
        rows.append([*cen, r])
    return np.array(rows)

def certify(support_path, queries, key):
    p = D / (key + '_queries.txt')
    with p.open('w') as h:
        h.write(str(len(queries)) + '\n')
        np.savetxt(h, queries, fmt='%.17g')
    z = subprocess.run([str(D / 'exact_clearance'), str(support_path), str(p)], capture_output=True, text=True, timeout=90)
    if z.returncode:
        raise ValueError('EXACT_CERT_EXIT ' + str(z.returncode) + ' ' + z.stderr[:100])
    out = json.loads(z.stdout)
    out.update(support_path=str(support_path), support_sha256=sha(support_path), query_path=str(p), query_sha256=sha(p), all_pass=all((q['pass'] and q.get('inside', True) for q in out['rows'])))
    dump(R / 'raw' / (key + '_EXACT.json'), out)
    return dict(path=str(R / 'raw' / (key + '_EXACT.json')), sha256=sha(R / 'raw' / (key + '_EXACT.json')), all_pass=out['all_pass'], balls=len(queries), scope='Exact rational triangle separation for outward-rounded covering radii; topological inside checked by exact CGAL when present')

def normals_exact(v, f, d):
    dq = [F(float(x)) for x in d]
    mn = None
    bad = 0
    for tri in v[f]:
        (a, b, c) = [[F(float(x)) for x in row] for row in tri]
        u = [b[i] - a[i] for i in range(3)]
        w = [c[i] - a[i] for i in range(3)]
        n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
        dot = sum((x * y for (x, y) in zip(n, dq)))
        bad += dot <= 0
        mn = dot if mn is None else min(mn, dot)
    return dict(pass_exact=bad == 0, bad_facets=bad, minimum_unnormalized_dot_exact=str(mn), direction=[str(x) for x in dq], scope='Exact stored-coordinate outward-normal witness on nonbasal facets. Monotone ring solid supplies isolated-tooth translation path; full arch remains separate.')
