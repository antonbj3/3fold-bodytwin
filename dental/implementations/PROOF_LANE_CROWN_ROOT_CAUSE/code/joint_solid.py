"""Joint distance construction; artificial closure is always labelled separately."""
from dental_release.paths import expand as _release_expand
from diagnose import *
import subprocess, datetime
from skimage.measure import marching_cubes
import igl
D = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/XREV-batch32/external/PROOF_LANE_CROWN_ROOT_CAUSE'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def meshwrite(p, v, f):
    with Path(p).open('w') as h:
        h.write(f'{len(v)} {len(f)}\n')
        np.savetxt(h, v, fmt='%.17g')
        np.savetxt(h, f, fmt='%d')

def meshread(p):
    with Path(p).open() as h:
        (nv, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=nv)
        f = np.loadtxt(h, dtype=int)
    return (v, f)

def close(rec):
    tri = load(rec['private_path'])['target']
    (v, inv) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
    f = inv.reshape(-1, 3)
    mm = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(mm, multibody=True)
    p = D / (rec['key'] + '_native.mesh')
    dest = D / (rec['key'] + '_closed.mesh')
    meshwrite(p, mm.vertices, mm.faces)
    proc = subprocess.run([str(D / 'close_native'), str(p), str(dest)], capture_output=True, text=True)
    if proc.returncode:
        raise ValueError('CGAL closure error ' + str(proc.returncode) + ' ' + proc.stderr[:200])
    info = json.loads(proc.stdout)
    (v, f) = meshread(dest)
    mm = trimesh.Trimesh(v, f, process=False)
    if mm.volume < 0:
        f = f[:, ::-1]
    source_set = {tuple(sorted((tuple(x) for x in t))) for t in tri}
    native = np.array([tuple(sorted((tuple(x) for x in t))) in source_set for t in v[f]])
    info.update(native_retained_count=int(native.sum()), source_count=len(tri), watertight=bool(mm.is_watertight), volume_mm3=float(abs(mm.volume)))
    dump(R / 'raw' / ('B_CLOSURE_' + rec['key'] + '.json'), info)
    if info['intersection_count'] or not info['closed'] or native.sum() != len(tri):
        raise ValueError('INVALID_SOURCE_SUPPORT ' + json.dumps(info))
    return (v, f, native, info)

def iso(field, lo, h):
    (v, f, _, _) = marching_cubes(field, level=0, spacing=(h,) * 3, allow_degenerate=False)
    v = v.astype(float) + lo
    m = trimesh.Trimesh(v, f, process=False)
    trimesh.repair.fix_normals(m, multibody=True)
    if m.volume < 0:
        m.invert()
    return (m.vertices, m.faces)

def generate(rec, closure_only=False):
    st = time.perf_counter()
    (v, f, native, info) = close(rec)
    if closure_only:
        return dict(key=rec['key'], status='VALID_CLOSURE', closure=info)
    pr = read(R / 'PREREG_B.json')['metrics']
    h = pr['grid_pitch_mm']
    r = pr['crown_level_mm']
    S = v[f[native]].copy()
    A = v[f[~native]].copy()
    lo = v.min(0) - 2 * h
    hi = v.max(0) + 2 * h
    axes = [np.arange(lo[k], hi[k] + h, h) for k in range(3)]
    shape = tuple(map(len, axes))
    size = int(np.prod(shape))
    if size > 16000000:
        raise ValueError('grid nodes exceed frozen budget: ' + str(size))
    sd = np.empty(size, np.float32)
    ds = sd.copy()
    dc = sd.copy()
    for start in range(0, size, 32768):
        j = np.arange(start, min(size, start + 32768))
        ijk = np.array(np.unravel_index(j, shape)).T
        pts = np.column_stack([axes[k][ijk[:, k]] for k in range(3)])
        (ss, _, _, _) = igl.signed_distance(pts, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        sd[j] = ss
        ds[j] = fast_nearest(S, pts)[1]
        dc[j] = fast_nearest(A, pts)[1]
    gap = 0.025 + 0.025 * np.clip(dc, 0, 1)
    crown = np.maximum(sd, ds - r).reshape(shape)
    prep = np.maximum(sd, r + gap - ds).reshape(shape)
    (cv, cf) = iso(crown, lo, h)
    (pv, pf) = iso(prep, lo, h)

    def roles_for(vv, ff, part):
        pts = vv[ff].mean(1)
        aa = fast_nearest(A, pts)[1]
        ss = fast_nearest(S, pts)[1]
        if part == 'crown':
            return np.where(ss < 0.5 * pr['crown_level_mm'], 0, np.where(aa < abs(ss - r), 2, 1))
        return np.where(aa < abs(ss - (r + 0.025 + 0.025 * np.clip(aa, 0, 1))), 2, 1)
    roles = roles_for(cv, cf, 'crown')
    proles = roles_for(pv, pf, 'prep')
    path = D / (rec['key'] + '_B.npz')
    np.savez_compressed(path, vertices=cv, faces=cf, roles=roles, prep_vertices=pv, prep_faces=pf, prep_roles=proles, native_triangles=S, closure_triangles=A)
    out = dict(key=rec['key'], family=rec['family'], status='GENERATED', mesh_path=str(path), mesh_sha256=sha(path), closure=info, grid_nodes=size, grid_shape=list(shape), seconds=time.perf_counter() - st)
    print('GENERATED', rec['key'], len(cv), len(cf), out['seconds'], flush=True)
    return out

def run():
    records = read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']
    keys = read(R / 'PREREG_B.json')['selection']['primary_keys']
    records = [r for r in records if '--all' in sys.argv or r['key'] in keys]
    rows = []
    st = time.perf_counter()
    for rec in records:
        try:
            out = generate(rec, '--closure-only' in sys.argv)
        except Exception as e:
            out = dict(key=rec['key'], family=rec['family'], status='REJECTED', reason=repr(e))
            print(out, flush=True)
        rows.append(out)
        dump(R / 'raw/B_GENERATION.json', rows)
        dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-root-cause', phase='B_CONSTRUCTING', latest_gate=out, next_operation='Freeze complete outputs before separate evaluation; retain closure failures', review_state='PENDING_INDEPENDENT_REVIEW'))
    if '--closure-only' not in sys.argv:
        p = R / 'FROZEN_PREDICTIONS_B.json'
        if p.exists():
            raise ValueError('Predictions already frozen')
        dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(R / 'PREREG_B.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, physical_predictions=None, scope='Digital geometry only; native exterior used as input, no held-out inference'))
        p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
if __name__ == '__main__':
    run()
