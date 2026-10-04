from fc_common import *
from scipy import ndimage
from skimage.measure import marching_cubes
import trimesh, resource
F = functional()

def classify(points, tri, margin, neighbors):
    sm = F.make_mesh(tri)
    ids = []
    for x in np.array_split(points, max(1, int(np.ceil(len(points) / 128)))):
        ids.extend(trimesh.proximity.closest_point(sm, x)[2])
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-12)
    nz = np.abs(normals[np.array(ids), 2])
    center = tri.reshape(-1, 3).mean(0)
    xy = points[:, :2] - center[:2]
    xy /= np.maximum(np.linalg.norm(xy, axis=1, keepdims=True), 1e-12)
    prox = np.zeros(len(points), bool)
    for n in neighbors:
        if not len(n):
            continue
        d = n.reshape(-1, 3).mean(0)[:2] - center[:2]
        d /= max(np.linalg.norm(d), 1e-12)
        prox |= xy @ d >= np.sqrt(0.5)
    labels = np.where(prox, 'proximal', 'buccolingual_axial')
    labels = np.where(nz >= 0.6, 'occlusal', labels)
    labels = np.where(points[:, 2] <= margin + 0.75, 'margin', labels)
    return labels

def run():
    st = time.perf_counter()
    records = read(V4 / 'payload/whole_inputs/RECORDS.json')
    res = read(V4 / 'results.json')
    rows = []
    oracles = []
    for md in records:
        key = md['key']
        a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
        p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
        margin = float(p['margin_z'])
        tri = a['source_triangles']
        tri = tri[tri.mean(1)[:, 2] >= margin]
        neigh = [a[s + '_triangles'] for s in ['mesial', 'distal']]
        for (track, folder) in [('R3', 'whole_scores'), ('R5', 'prep_scores')]:
            for old in res['rounds'][track]['rows']:
                if old['key'] != key or old.get('full_surface_p95_mm') is None:
                    continue
                t = npz(V4 / 'payload' / folder / old['participant'] / (key + '.npz'))
                for (direction, pp, dd) in [('generated_to_source', t['exterior_samples'], t['exterior_to_source_mm']), ('source_to_generated', t['source_samples'], t['source_to_exterior_mm'])]:
                    lab = classify(pp, tri, margin, neigh)
                    p95 = np.quantile(dd, 0.95)
                    tail = dd >= p95
                    for name in ['margin', 'occlusal', 'proximal', 'buccolingual_axial']:
                        mask = lab == name
                        rows.append(dict(key=key, family=md['family'], participant=old['participant'], track=track, direction=direction, region=name, n=int(mask.sum()), tail_count=int(np.sum(mask & tail)), tail_total=int(tail.sum()), p95_mm=float(np.quantile(dd[mask], 0.95)) if mask.any() else None, global_p95_mm=float(p95), resolution='PER_SURFACE_REGION'))
        h = float(p['step'])
        zz = p['origin'][2] + h * np.arange(p['preparation'].shape[2])
        cut = np.broadcast_to(margin - zz[None, None, :], p['preparation'].shape)
        outer = a['phi_source']
        field = np.maximum.reduce([outer, -p['cavity'], cut])
        (v, f, _, _) = marching_cubes(field, 0.001, spacing=(h,) * 3, allow_degenerate=False)
        v += p['origin']
        pts = v[f].mean(1)
        ix = ((pts - p['origin']) / h).T
        roles = np.argmax(np.stack([ndimage.map_coordinates(outer, ix, order=1), ndimage.map_coordinates(-p['cavity'], ix, order=1), margin - pts[:, 2]]), axis=0).astype(np.int8)
        if trimesh.Trimesh(v, f, process=False).volume < 0:
            f = f[:, ::-1]
        d = DATA / 'diagnostic_oracles' / key
        d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
        rr = F.score_mesh(dict(md, participant='known_source_insertion_envelope_ORACLE'), v, f, roles)
        oracles.append(rr)
        dump(ROOT / 'raw/R1_ORACLES_CHECKPOINT.json', oracles)
        print(key, 'oracle p95', round(rr['reconstruction_p95_mm'], 4), flush=True)
    out = dict(round='R1', claim_type='capability', external_referent=REFERENT, oracle_pass=sum((r['reconstruction_p95_mm'] <= 0.35 for r in oracles)), oracle_requested=len(records), oracle_rows=oracles, regional_rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'rounds/R1.json', out)
    state('R1_DECIDED', f"known-source envelope passes {out['oracle_pass']}/{len(records)}", 'Construct full oriented-surface prior conditional on public preparation; freeze before score')
if __name__ == '__main__':
    run()
