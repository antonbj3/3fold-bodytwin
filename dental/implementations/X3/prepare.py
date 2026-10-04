from crownbench import *
import argparse, resource

def train():
    start = time.perf_counter()
    split = json.load(open(P / 'SPLIT.json'))
    acc = {}
    metadata = {}
    raw = []
    state('TRAIN_PRIOR', 'PREREG_R1_FROZEN', 'finish train-only signed fields')
    for rec in split['scans']:
        if rec['split'] != 'train':
            continue
        (v, f) = load_obj(rec['mesh'])
        j = json.load(open(rec['labels']))
        labels = np.asarray(j['labels'], dtype=np.int16)
        assert len(v) == len(labels), 'VERTEX_ORDER_OR_LABEL_LENGTH'
        teeth = tooth_dict(v, f, labels)
        n = axis(teeth)
        for i in TARGETS:
            fields = []
            border = []
            for q in [1, 2] if rec['jaw'] == 'upper' else [3, 4]:
                k = q * 10 + i
                if k not in teeth:
                    continue
                (tv, tf) = teeth[k]
                if len(tf) < 200:
                    continue
                try:
                    F = local_frame(teeth, k, n)
                    lv = (tv - tv.mean(0)) @ F
                    if q in [2, 4]:
                        lv[:, 1] *= -1
                        tf = tf[:, ::-1]
                    phi = signed_grid(lv, tf)
                    (_, ids) = boundary(lv, tf)
                    fields.append(phi)
                    border.append(float(np.median(lv[ids, 2])))
                    raw.append({'case': rec['case'], 'jaw': rec['jaw'], 'fdi': k, 'vertices': len(tv), 'faces': len(tf), 'cap_boundary_vertices': len(ids)})
                    if not (DATA / 'train_example.npz').exists():
                        write_surface(DATA / 'train_example.npz', tv, teeth[k][1])
                        put(DATA / 'train_example_frame.json', {'F': F.tolist(), 'center': tv.mean(0).tolist()})
                except Exception as e:
                    raw.append({'case': rec['case'], 'jaw': rec['jaw'], 'fdi': k, 'failure': str(e)})
            if fields:
                key = rec['jaw'] + '_' + str(i)
                field = np.mean(fields, axis=0)
                if key not in acc:
                    acc[key] = np.zeros(SHAPE, dtype=np.float64)
                    metadata[key] = {'case_count': 0, 'boundary_z': []}
                acc[key] += field
                metadata[key]['case_count'] += 1
                metadata[key]['boundary_z'].extend(border)
        print('TRAIN', rec['case'], rec['jaw'], flush=True)
    for (key, a) in acc.items():
        phi = (a / metadata[key]['case_count']).astype(np.float32)
        (pv, pf) = zero_mesh(phi)
        tri = pv[pf]
        normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        normal /= np.maximum(np.linalg.norm(normal, axis=1, keepdims=True), 1e-12)
        bz = float(np.median(metadata[key]['boundary_z']))
        keep = ~((tri.mean(1)[:, 2] < bz + 0.2) & (normal[:, 2] < -0.35))
        pf = pf[keep]
        write_surface(DATA / 'model' / f'{key}.npz', pv, pf)
        np.savez_compressed(DATA / 'model' / f'{key}_field.npz', phi=phi, origin=LO, h=H)
        metadata[key]['boundary_median_z'] = bz
        metadata[key]['cap_triangles_removed'] = int((~keep).sum())
        del metadata[key]['boundary_z']
    (v, f) = read_surface(DATA / 'train_example.npz')
    e = json.load(open(DATA / 'train_example_frame.json'))
    (sv, sf, phi, detail) = sdf_roundtrip(v, f, np.array(e['F']), np.array(e['center']))
    pts = sampled(v, f, 2048, 6103)
    ds = distances(pts, sv, sf)
    n = np.array(e['F'])[:, 2]
    err = float(np.quantile(ds[pts @ n >= np.quantile(pts @ n, 0.65)], 0.95))
    detail.update(train_topband_p95_mm=err, passes=err <= 0.3)
    put(P / 'TRAIN_FIT.json', {'metadata': metadata, 'samples': raw, 'seconds': time.perf_counter() - start, 'maxrss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'roundtrip': detail, 'model_manifest': [{'path': str(p), 'sha256': sha(p)} for p in (DATA / 'model').glob('*.npz')]})
    print('TRAIN_DONE', time.perf_counter() - start, 'roundtrip', err, flush=True)
    state('TRAIN_COMPLETE', {'roundtrip': detail, 'data_bytes': disk_guard()}, 'prepare masked public R1 contexts')

def contexts(group):
    start = time.perf_counter()
    rows = []
    for rec in json.load(open(P / 'SPLIT.json'))['scans']:
        if rec['split'] != group:
            continue
        assert sha(pathlib.Path(rec['mesh'])) == rec['mesh_sha256']
        assert sha(pathlib.Path(rec['labels'])) == rec['labels_sha256']
        (v, f) = load_obj(rec['mesh'])
        labels = np.asarray(json.load(open(rec['labels']))['labels'], dtype=np.int16)
        assert len(v) == len(labels)
        teeth = tooth_dict(v, f, labels)
        q = 1 if rec['jaw'] == 'upper' else 3
        for i in TARGETS:
            k = q * 10 + i
            key = rec['case'] + '_' + rec['jaw'] + '_' + str(k)
            row = {'key': key, 'case': rec['case'], 'jaw': rec['jaw'], 'fdi': k, 'split': group, 'eligible': False, 'antagonist_contact': 'UNKNOWN_NO_VERIFIED_BITE_TRANSFORM'}
            req = [k, swap(k), k - 1, k + 1]
            if any((z not in teeth or len(teeth[z][1]) < 200 for z in req)):
                row['reason'] = 'MISSING_TARGET_HOMOLOG_OR_NEIGHBOR_OR_FACE_COUNT'
            else:
                (mv, mf, ml) = mask(v, f, labels, k)
                vis = tooth_dict(mv, mf, ml)
                try:
                    reflection(vis, q)
                    assert k not in np.unique(ml), 'TARGET_LEAKAGE'
                    context = DATA / 'public' / group / f'{key}.npz'
                    context.parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(context, vertices=mv, faces=mf, labels=ml, target_fdi=k)
                    hidden = DATA / 'hidden' / group / f'{key}.npz'
                    write_surface(hidden, *teeth[k])
                    row.update(eligible=True, context=str(context), context_sha256=sha(context), reference=str(hidden), reference_sha256=sha(hidden), visible_vertices=len(mv), visible_faces=len(mf), removed_vertices=int((labels == k).sum()), removed_incident_faces=len(f) - len(mf))
                except Exception as e:
                    row['reason'] = str(e)
            rows.append(row)
        print('MASK', group, rec['case'], rec['jaw'], flush=True)
    public = [{k: v for (k, v) in r.items() if k not in ('reference', 'reference_sha256', 'removed_vertices', 'removed_incident_faces')} for r in rows]
    put(DATA / 'public' / f'{group}_INDEX.json', public)
    put(DATA / 'hidden' / f'{group}_INDEX.json', rows)
    put(P / f'PREPARATION_{group}.json', {'scheduled': len(rows), 'eligible': sum((r['eligible'] for r in rows)), 'seconds': time.perf_counter() - start, 'data_bytes': disk_guard(), 'maxrss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    state(group + '_CONTEXTS_READY', {'scheduled': len(rows), 'eligible': sum((r['eligible'] for r in rows))}, 'generate predictions without loading hidden files')
if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('action', choices=['train', 'test_R1', 'reserve_R2'])
    args = a.parse_args()
    if args.action == 'train':
        train()
    else:
        contexts(args.action)
