from construct_g import *

def split_chosen(v, f, chosen):
    verts = v.tolist()
    mid = {}
    for (a, b) in sorted(chosen):
        mid[a, b] = len(verts)
        verts.append(((v[a] + v[b]) / 2).tolist())
    out = []
    for face in f:
        mids = [mid.get(tuple(sorted((int(face[j]), int(face[(j + 1) % 3]))))) for j in range(3)]
        n = sum((x is not None for x in mids))
        if n == 0:
            out.append(face.tolist())
        elif n == 1:
            j = next((j for j in range(3) if mids[j] is not None))
            (a, b, c) = np.roll(face, -j)
            m = mids[j]
            out.extend([[a, m, c], [m, b, c]])
        elif n == 2:
            j = next((j for j in range(3) if mids[j] is not None and mids[(j - 1) % 3] is not None))
            a = face[j]
            b = face[(j + 1) % 3]
            c = face[(j - 1) % 3]
            u = mids[j]
            w = mids[(j - 1) % 3]
            out.extend([[a, u, w], [w, u, b], [w, b, c]])
        else:
            (a, b, c) = face
            (u, w, z) = mids
            out.extend([[a, u, z], [u, b, w], [z, w, c], [u, w, z]])
    return (np.array(verts), np.array(out, int))

def make(m, key):
    (pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
    rim = pv[loops(trimesh.Trimesh(pv, pf, process=False))[0]]
    normal = trimesh.Trimesh(pv, pf, process=False).face_normals
    path = D6 / (key + '_G_closure_1.5_after.mesh')
    if not path.exists():
        path = D6 / (key + '_I_closure_1.5_after.mesh')
    with path.open() as h:
        (nv, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=nv)
        f = np.loadtxt(h, dtype=int)
    source_v = {tuple(p) for p in pv}
    source_t = {tuple(sorted((tuple(p) for p in t))) for t in pv[pf]}
    history = []
    for iteration in range(11):
        si = intersect(v, f, key + '_J_' + str(iteration))
        keep = retained_faces(v, f, pv, pf)
        cm = trimesh.Trimesh(v, f, process=False)
        history.append(dict(iteration=iteration, intersections=si.get('count'), source_retained=int(keep.sum()), watertight=cm.is_watertight, vertices=len(v), faces=len(f)))
        save(R6 / 'raw' / ('J_CLOSURE_' + key + '.json'), history)
        if not keep.all() or not cm.is_watertight:
            raise ValueError('J lost immutable source or closure')
        if si['status'] == 'PASS':
            isorig = np.array([tuple(sorted((tuple(p) for p in t))) in source_t for t in v[f]])
            dest = D6 / (key + '_J_closure.npz')
            np.savez_compressed(dest, vertices=v, faces=f)
            return (v, f, pv, pf, rim, isorig, history)
        if iteration == 10:
            break
        isorig = np.array([tuple(sorted((tuple(p) for p in t))) in source_t for t in v[f]])
        sourceedges = {tuple(sorted(e)) for face in f[isorig] for e in zip(face, np.roll(face, -1))}
        bad = set((int(i) for p in si['pairs'] for i in p if not isorig[i]))
        chosen = set()
        for i in bad:
            ee = [tuple(sorted(e)) for e in zip(f[i], np.roll(f[i], -1))]
            allowed = [e for e in ee if e not in sourceedges]
            if allowed:
                chosen.add(max(allowed, key=lambda e: np.linalg.norm(v[e[1]] - v[e[0]])))
        oldnv = len(v)
        affected = set((int(x) for i in bad for x in f[i]))
        (v, f) = split_chosen(v, f, chosen)
        affected.update(range(oldnv, len(v)))
        free = np.array([i for i in sorted(affected) if tuple(v[i]) not in source_v], int)
        if not len(free):
            raise ValueError('J no mutable collision degrees')
        (q, d, j) = fast_nearest(pv[pf], v[free])
        sign = np.einsum('ij,ij->i', v[free] - q, normal[j])
        active = (d < 0.25) & (sign > -0.1)
        v[free[active]] = q[active] - 0.1 * normal[j[active]]
    raise ValueError('J collision budget exhausted ' + str(history[-1]))
if __name__ == '__main__':
    (a, b) = first()
    st = time.perf_counter()
    try:
        (v, f, pv, pf, rim, isorig, hi) = make(npz(a['mesh_path']), a['key'])
        out = dict(status='VALID_CLOSURE', history=hi)
    except Exception as e:
        out = dict(status='REJECTED', reason=repr(e))
    out.update(seconds=time.perf_counter() - st, claim_type='capability')
    save(R6 / 'RESULTS_J.json', out)
    print(clean(out))
