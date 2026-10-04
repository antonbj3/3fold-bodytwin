"""Fixed RAS pose, projected surface contact and conditional tangent region load.
Raster operation adapted from predecessor P/code/occlusion.py (raster function).
No virtual jaw transform. Physical force and true Bits2Bites labels remain unknown.
"""
from dental_release.paths import expand as _release_expand
import hashlib, json, time, resource, zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, nnls
H = Path(__file__).resolve().parents[1]
DATA = Path(__import__('os').environ.get('X2_DATA_DIR', _release_expand('@DENTAL_WORK_ROOT@/X2_occlusion_b2b')))
ZIP = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'))
TDS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/Teeth3DS'))
REGIONS = ['I_R', 'C_R', 'P_R', 'M_R', 'I_L', 'C_L', 'P_L', 'M_L']
TYPE_TO_REGION = np.array([0, 0, 0, 1, 2, 2, 3, 3, 3])

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False, default=lambda a: a.item() if isinstance(a, np.generic) else a.tolist() if isinstance(a, np.ndarray) else str(a)) + '\n')
    tmp.replace(p)

def state(stage, gate, next_op, **kwargs):
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X2-occlusion-b2b', status=stage, latest_gate=gate, next_operation=next_op, updated_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), **kwargs))

def features(v, force_ras=False):
    (lo, hi) = np.quantile(v[:, :2], [0.01, 0.99], axis=0)
    mid = (lo + hi) / 2
    scale = np.maximum((hi - lo) / 2, 1e-09)
    xx = (v[:, 0] - mid[0]) / scale[0]
    center = v[np.abs(xx) < 0.25, 1]
    side = v[np.abs(xx) > 0.65, 1]
    ysign = 1 if force_ras or np.median(center) >= np.median(side) else -1
    y = (v[:, 1] - mid[1]) * ysign / scale[1]
    return (np.c_[np.abs(xx), y], dict(x_center_mm=float(mid[0]), y_center_mm=float(mid[1]), xy_halfspan_mm=scale.tolist(), anterior_sign=int(ysign)))

def predict_regions(v, model, ras=False):
    (X, frame) = features(v, ras)
    proto = np.asarray(model['prototypes'])
    pred = np.sum((X[:, None, :] - proto[None, :, :]) ** 2, axis=2).argmin(1) + 1
    r = TYPE_TO_REGION[pred] + 4 * (v[:, 0] < frame['x_center_mm'])
    return (r, pred, frame)

def load_obj(p):
    vv = []
    with open(p) as f:
        for l in f:
            if l.startswith('v '):
                vv.append([float(x) for x in l.split()[1:4]])
    return np.asarray(vv)

def teeth_files(jaw, split):
    names = (TDS / 'Teeth3DS_train_test_split' / f'{split}_{jaw}.txt').read_text().splitlines()
    out = []
    for n in names:
        n = n.strip()
        pid = n.rsplit('_', 1)[0]
        for part in range(1, 8):
            p = TDS / f'data_part_{part}' / jaw / pid / (n + '.obj')
            if p.exists() and p.with_suffix('.json').exists():
                out.append((pid, p))
                break
    return sorted(out, key=lambda x: hashlib.sha256(x[0].encode()).hexdigest())

def train_and_validate():
    start = time.perf_counter()
    models = {}
    rows = []
    cost = []
    for jaw in ['upper', 'lower']:
        fit = []
        paths = teeth_files(jaw, 'training')[:50]
        for (pid, p) in paths:
            v = load_obj(p)
            lab = np.asarray(json.loads(p.with_suffix('.json').read_text())['labels'])
            (X, frame) = features(v)
            if len(lab) != len(v):
                raise ValueError('Labels mismatch ' + pid)
            for t in range(1, 9):
                for quadrant in [1, 2] if jaw == 'upper' else [3, 4]:
                    ids = lab == quadrant * 10 + t
                    if ids.sum() >= 50:
                        fit.append((t, X[ids].mean(0)))
            cost.append(dict(patient=pid, jaw=jaw, role='fit', mesh_path=str(p), mesh_sha256=sha(p), labels_sha256=sha(p.with_suffix('.json'))))
        model = dict(prototypes=[np.mean([x for (t0, x) in fit if t0 == t], axis=0).tolist() for t in range(1, 9)], training_arches=len(paths))
        models[jaw] = model
        confusion = np.zeros((4, 4), int)
        for (pid, p) in teeth_files(jaw, 'testing')[:30]:
            v = load_obj(p)
            lab = np.asarray(json.loads(p.with_suffix('.json').read_text())['labels'])
            (r, ty, frame) = predict_regions(v, model)
            ok = lab > 0
            true = TYPE_TO_REGION[lab[ok] % 10]
            pred = r[ok] % 4
            cf = np.bincount(true * 4 + pred, minlength=16).reshape(4, 4)
            confusion += cf
            rows.append(dict(patient=pid, jaw=jaw, confusion=cf.tolist(), accuracy=float(np.trace(cf) / cf.sum()), frame=frame, mesh_path=str(p), mesh_sha256=sha(p), labels_sha256=sha(p.with_suffix('.json'))))
        acc = np.diag(confusion) / np.maximum(confusion.sum(1), 1)
        ap_error = (confusion[:2, 2:].sum() + confusion[2:, :2].sum()) / confusion.sum()
        model['validation'] = dict(confusion=confusion.tolist(), per_region_recall=acc.tolist(), macro_recall=float(acc.mean()), anterior_posterior_error=float(ap_error), n_arches=sum((x['jaw'] == jaw for x in rows)), gate=bool(acc.mean() >= 0.85 and ap_error <= 0.05))
        print('region validation', jaw, model['validation'], flush=True)
    dump(H / 'raw/segmentation_model.json', models)
    dump(H / 'raw/segmentation_validation.json', dict(rows=rows, training_inputs=cost, wall_s=time.perf_counter() - start, models=models))
    return models

def read_stl(z, name):
    raw = z.read(name)
    n = int(np.frombuffer(raw[80:84], dtype='<u4')[0])
    dtype = np.dtype([('normal', '<f4', (3,)), ('triangle', '<f4', (3, 3)), ('attr', '<u2')])
    if len(raw) != 84 + 50 * n:
        raise ValueError('Non-binary STL length mismatch')
    t = np.frombuffer(raw, offset=84, count=n, dtype=dtype)['triangle'].astype(float)
    return (t, hashlib.sha256(raw).hexdigest())

def raster(tri, h, origin, shape, upper):
    n = len(tri)
    ij = np.floor((tri[:, :, :2].min(1) - origin) / h - 0.5).astype(int)
    stop = np.ceil((tri[:, :, :2].max(1) - origin) / h - 0.5).astype(int)
    spans = stop - ij
    grid = np.full(np.prod(shape), np.inf if upper else -np.inf)
    owner = np.full(len(grid), -1, int)

    def update(ids, loc, xy):
        t = tri[ids]
        (a, b, c) = (t[:, 0, :2], t[:, 1, :2], t[:, 2, :2])
        den = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
        ok = np.abs(den) > 1e-12
        w1 = np.divide((b[:, 1] - c[:, 1]) * (xy[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (xy[:, 1] - c[:, 1]), den, out=np.zeros_like(den), where=ok)
        w2 = np.divide((c[:, 1] - a[:, 1]) * (xy[:, 0] - c[:, 0]) + (a[:, 0] - c[:, 0]) * (xy[:, 1] - c[:, 1]), den, out=np.zeros_like(den), where=ok)
        w3 = 1 - w1 - w2
        mask = ok & (w1 >= -1e-08) & (w2 >= -1e-08) & (w3 >= -1e-08) & (loc[:, 0] >= 0) & (loc[:, 1] >= 0) & (loc[:, 0] < shape[0]) & (loc[:, 1] < shape[1])
        ids = ids[mask]
        if not len(ids):
            return
        ix = loc[mask, 0] * shape[1] + loc[mask, 1]
        height = (w1 * t[:, 0, 2] + w2 * t[:, 1, 2] + w3 * t[:, 2, 2])[mask]
        order = np.lexsort((ids, height if upper else -height, ix))
        ordered = ix[order]
        keep = np.r_[True, np.diff(ordered) != 0]
        ix = ordered[keep]
        selected = order[keep]
        height = height[selected]
        ids = ids[selected]
        better = height < grid[ix] if upper else height > grid[ix]
        grid[ix[better]] = height[better]
        owner[ix[better]] = ids[better]
    largest = spans.max(1)
    previous = -1
    for cap in [2, 4, 8, 16]:
        group = np.flatnonzero((largest > previous) & (largest <= cap))
        previous = cap
        if not len(group):
            continue
        span = spans[group].max(0)
        for dx in range(int(span[0]) + 1):
            for dy in range(int(span[1]) + 1):
                ids = group[(spans[group, 0] >= dx) & (spans[group, 1] >= dy)]
                loc = ij[ids] + [dx, dy]
                update(ids, loc, origin + (loc + 0.5) * h)
    for i in np.flatnonzero(largest > 16):
        (xx, yy) = np.meshgrid(np.arange(ij[i, 0], stop[i, 0] + 1), np.arange(ij[i, 1], stop[i, 1] + 1), indexing='ij')
        loc = np.c_[xx.ravel(), yy.ravel()]
        update(np.full(len(loc), i), loc, origin + (loc + 0.5) * h)
    return (grid, owner)

def make_map(U, L, h):
    lo = np.minimum(U[:, :, :2].min((0, 1)), L[:, :, :2].min((0, 1)))
    hi = np.maximum(U[:, :, :2].max((0, 1)), L[:, :, :2].max((0, 1)))
    origin = np.floor(lo / h) * h - 2 * h
    shape = np.ceil((hi - origin) / h).astype(int) + 3
    (uz, ui) = raster(U, h, origin, shape, True)
    (lz, li) = raster(L, h, origin, shape, False)
    ok = np.isfinite(uz) & np.isfinite(lz)
    ids = np.flatnonzero(ok)
    xy = origin + (np.c_[ids // shape[1], ids % shape[1]] + 0.5) * h
    return dict(gap=uz[ok] - lz[ok], upper_face=ui[ok], lower_face=li[ok], xy=xy, upper_z=uz[ok], lower_z=lz[ok], origin=origin, shape=shape)

def tangent(pairs, k):
    C = np.zeros((len(pairs), 16))
    C[np.arange(len(pairs)), pairs[:, 0]] = 1
    C[np.arange(len(pairs)), pairs[:, 1] + 8] = 1
    if len(pairs) == 0:
        return None
    kr = k / k.mean()
    r = minimize(lambda u: 0.5 * np.dot(kr * u, u), np.ones(16) * 0.5, jac=lambda u: kr * u, method='SLSQP', bounds=[(0, None)] * 16, constraints=[dict(type='ineq', fun=lambda u: C @ u - 1, jac=lambda u: C)], options=dict(ftol=1e-13, maxiter=250))
    u = r.x
    F = k * u
    slack = C @ u - 1
    active = slack < 1e-07
    lam = np.zeros(len(C))
    lam[active] = nnls(C[active].T, kr * u, maxiter=1000)[0]
    err = float(np.max(np.abs(C.T @ lam - kr * u)))
    comp = float(np.max(np.abs(lam * slack)))
    conservation = float(abs(F[:8].sum() - F[8:].sum()) / max(F[:8].sum(), 1e-15))
    w = np.r_[F[:8] / F[:8].sum(), F[8:] / F[8:].sum()] * 100
    return dict(shares_pp=w.tolist(), unit_increment_displacements=u.tolist(), conditional_force_per_unit_increment_N=F.tolist(), pairs=pairs.tolist(), valid=bool(slack.min() >= -1e-06 and err <= 1e-06 and (comp <= 1e-06) and (conservation <= 1e-07)), kkt_residual=err, complementarity=comp, conservation_relative=conservation, minimum_slack=float(slack.min()))

def analyze_map(m, Uregions, Lregions):
    gap = m['gap']
    ur = Uregions[m['upper_face']]
    lr = Lregions[m['lower_face']]
    rows = []
    k = np.tile([2, 1, 2, 2, 2, 1, 2, 2], 2) * np.sqrt(500 * 1130)
    for band in [0.05, 0.1, 0.2]:
        ok = gap <= band
        up = np.bincount(ur[ok], minlength=8)
        lo = np.bincount(lr[ok], minlength=8)
        if ok.sum() == 0:
            rows.append(dict(band_mm=band, n_pixels=0, mechanics=None, area_shares_pp=None))
            continue
        pairs = np.unique(np.c_[ur[ok], lr[ok]], axis=0)
        start = time.perf_counter()
        mech = tangent(pairs, k)
        ms = time.perf_counter() - start
        start = time.perf_counter()
        area = np.r_[up, lo] / ok.sum() * 100
        ars = time.perf_counter() - start
        rows.append(dict(band_mm=band, n_pixels=int(ok.sum()), area_shares_pp=area.tolist(), mechanics=mech, mechanics_query_s=ms, area_query_s=ars, projected_area_mm2=float(ok.sum() * m['grid_mm'] ** 2), negative_gap_pixels=int((gap[ok] < 0).sum())))
    return rows

def run_cases(models, max_cases=200, refinement=False):
    import csv, io
    DATA.mkdir(parents=True, exist_ok=True)
    if sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())) > 3000000000:
        raise RuntimeError('Lane disk limit')
    rows = []
    z = zipfile.ZipFile(ZIP)
    annotations = {int(r['Patient']): r for r in csv.DictReader(io.StringIO(z.read('Bits2Bites/data/Annotations.csv').decode('utf-8-sig')))}
    names = {int(n.split('/')[-2]): n for n in z.namelist() if n.endswith('upper.stl')}
    for case in range(1, max_cases + 1):
        outfile = H / 'raw/cases' / f'{case:03d}.json'
        if outfile.exists():
            rows.append(json.loads(outfile.read_text()))
            continue
        state('R1_GEOMETRY_RUNNING', 'Region transfer assessed; measured-force identification still UNKNOWN', f'Fixed-pose raster and contact graph case{case}', completed_cases=len(rows))
        st = time.perf_counter()
        cpu = time.process_time()
        (U, uh) = read_stl(z, names[case])
        (L, lh) = read_stl(z, names[case].replace('upper.stl', 'lower.stl'))
        decode_s = time.perf_counter() - st
        st = time.perf_counter()
        (ur, ut, uf) = predict_regions(U.mean(1), models['upper'], ras=True)
        (lr, lt, lf) = predict_regions(L.mean(1), models['lower'], ras=True)
        seg_s = time.perf_counter() - st
        st = time.perf_counter()
        m = make_map(U, L, 0.2)
        m['grid_mm'] = 0.2
        map_s = time.perf_counter() - st
        metrics = analyze_map(m, ur, lr)
        row = dict(case=case, annotation=annotations[case], upper_member=names[case], upper_sha256=uh, lower_sha256=lh, upper_frame=uf, lower_frame=lf, upper_triangles=len(U), lower_triangles=len(L), minimum_sampled_vertical_gap_mm=float(m['gap'].min()), gap_quantiles_mm=np.quantile(m['gap'], [0, 0.01, 0.05, 0.5]).tolist(), metrics=metrics, decode_s=decode_s, segmentation_s=seg_s, map_s=map_s, wall_s=time.perf_counter() - st + decode_s + seg_s, cpu_s=time.process_time() - cpu, rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
        np.savez_compressed(DATA / f'{case:03d}_map.npz', **{key: val for (key, val) in m.items() if key != 'grid_mm'}, upper_region=ur[m['upper_face']], lower_region=lr[m['lower_face']])
        row['map_path'] = str(DATA / f'{case:03d}_map.npz')
        row['map_sha256'] = sha(row['map_path'])
        dump(outfile, row)
        rows.append(row)
        print('case', case, 'gap', round(row['minimum_sampled_vertical_gap_mm'], 4), 'seconds', round(row['wall_s'], 2), flush=True)
    return rows
if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument('stage', choices=['fit', 'cases'])
    a.add_argument('--max-cases', type=int, default=200)
    args = a.parse_args()
    if args.stage == 'fit':
        train_and_validate()
    else:
        run_cases(json.loads((H / 'raw/segmentation_model.json').read_text()), args.max_cases)
